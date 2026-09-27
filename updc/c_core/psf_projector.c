#define _POSIX_C_SOURCE 200809L
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "overlap_geometry.h"

typedef struct {
    char input[1024];
    char vertices[1024];
    char output[1024];
    int width;
    int height;
    int out_width;
    int out_height;
    double input_sum;
    double output_sum;
    double flux_rel_error;
    uint64_t overlap_count;
} Job;

static Job *jobs = NULL;
static int njobs = 0;

static void die(const char *message) {
    fprintf(stderr, "UPDC_PSF_PROJECT: %s\n", message);
    exit(EXIT_FAILURE);
}

static void *xmalloc(size_t bytes) {
    void *p = malloc(bytes);
    if (!p) die("out of memory");
    return p;
}

static void *xcalloc(size_t count, size_t bytes) {
    void *p = calloc(count, bytes);
    if (!p) die("out of memory");
    return p;
}

static void read_f32_exact(const char *path, float *out, size_t count) {
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open PSF raw input");
    if (fread(out, sizeof(float), count, f) != count) die("short PSF raw input");
    if (fgetc(f) != EOF) die("PSF raw input has trailing bytes");
    fclose(f);
}

static void read_f64_exact(const char *path, double *out, size_t count) {
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open PSF vertex map");
    if (fread(out, sizeof(double), count, f) != count) die("short PSF vertex map");
    if (fgetc(f) != EOF) die("PSF vertex map has trailing bytes");
    fclose(f);
}

static double sum_f32(const float *data, size_t count) {
    double total = 0.0;
    for (size_t i = 0; i < count; ++i) {
        if (!isfinite(data[i]) || data[i] < 0.0f) die("PSF contains negative/non-finite pixels");
        total += (double)data[i];
    }
    return total;
}

static int symmetric_half_extent(double max_abs) {
    /* A target pixel with integer center t covers [t-0.5,t+0.5].
     * The small tolerance prevents roundoff at exact half-integer boundaries
     * from forcing an unnecessary extra one-pixel ring. */
    if (!(max_abs >= 0.0) || !isfinite(max_abs)) die("invalid PSF mapped extent");
    double needed = max_abs - 0.5 - 1.0e-7;
    if (needed <= 0.0) return 0;
    if (needed > 1000000.0) die("unreasonable PSF mapped extent");
    return (int)ceil(needed);
}

static void project_one(Job *job, int index) {
    const int iw = job->width, ih = job->height;
    const size_t nin = (size_t)iw * ih;
    const size_t nv = 2 * (size_t)(iw + 1) * (ih + 1);
    float *input = xmalloc(nin * sizeof(float));
    double *vertices = xmalloc(nv * sizeof(double));
    read_f32_exact(job->input, input, nin);
    read_f64_exact(job->vertices, vertices, nv);
    job->input_sum = sum_f32(input, nin);
    if (!(job->input_sum > 0.0)) die("PSF input sum is zero");

    double max_abs_x = 0.0, max_abs_y = 0.0;
    const size_t nvertex = (size_t)(iw + 1) * (ih + 1);
    for (size_t i = 0; i < nvertex; ++i) {
        const double x = vertices[2 * i], y = vertices[2 * i + 1];
        if (!isfinite(x) || !isfinite(y)) die("non-finite PSF mapped vertex");
        const double ax = fabs(x), ay = fabs(y);
        if (ax > max_abs_x) max_abs_x = ax;
        if (ay > max_abs_y) max_abs_y = ay;
    }
    const int hx = symmetric_half_extent(max_abs_x);
    const int hy = symmetric_half_extent(max_abs_y);
    const int ow = 2 * hx + 1, oh = 2 * hy + 1;
    if (ow < 1 || oh < 1) die("invalid PSF output shape");
    const size_t nout = (size_t)ow * oh;
    float *output = xcalloc(nout, sizeof(float));

    uint64_t overlaps = 0;
    UpdcPoint quad[4];
    for (int y = 0; y < ih; ++y) for (int x = 0; x < iw; ++x) {
        const size_t p = (size_t)y * iw + x;
        const float value = input[p];
        if (value == 0.0f) continue;
        const size_t a = 2 * ((size_t)y * (iw + 1) + x);
        const size_t b = a + 2 * (iw + 1);
        quad[0] = (UpdcPoint){vertices[a], vertices[a + 1]};
        quad[1] = (UpdcPoint){vertices[b], vertices[b + 1]};
        quad[2] = (UpdcPoint){vertices[b + 2], vertices[b + 3]};
        quad[3] = (UpdcPoint){vertices[a + 2], vertices[a + 3]};
        const float full_area = updc_polygon_area(quad, 4);
        if (!isfinite(full_area) || !(full_area > 1.0e-12f)) continue;

        double minx = quad[0].x, maxx = quad[0].x;
        double miny = quad[0].y, maxy = quad[0].y;
        for (int j = 1; j < 4; ++j) {
            if (quad[j].x < minx) minx = quad[j].x;
            if (quad[j].x > maxx) maxx = quad[j].x;
            if (quad[j].y < miny) miny = quad[j].y;
            if (quad[j].y > maxy) maxy = quad[j].y;
        }
        int x0 = (int)floor(minx + 0.5), x1 = (int)floor(maxx + 0.5);
        int y0 = (int)floor(miny + 0.5), y1 = (int)floor(maxy + 0.5);
        if (x0 < -hx) x0 = -hx;
        if (x1 >  hx) x1 =  hx;
        if (y0 < -hy) y0 = -hy;
        if (y1 >  hy) y1 =  hy;
        for (int ty = y0; ty <= y1; ++ty) for (int tx = x0; tx <= x1; ++tx) {
            const float area = updc_overlap_area(quad, tx, ty);
            if (!(area > 1.0e-8f)) continue;
            const float frac = area / full_area;
            const size_t q = (size_t)(ty + hy) * ow + (size_t)(tx + hx);
            output[q] = fmaf(value, frac, output[q]);
            ++overlaps;
        }
    }

    job->output_sum = sum_f32(output, nout);
    job->flux_rel_error = (job->output_sum - job->input_sum) / job->input_sum;
    job->overlap_count = overlaps;
    job->out_width = ow;
    job->out_height = oh;
    if (fabs(job->flux_rel_error) > 5.0e-5) {
        fprintf(stderr,
                "UPDC_PSF_PROJECT: exposure %d flux mismatch input=%.17g output=%.17g rel=%.6g\n",
                index + 1, job->input_sum, job->output_sum, job->flux_rel_error);
        die("PSF polygon projection is not flux-conserving");
    }

    FILE *f = fopen(job->output, "wb");
    if (!f) die("cannot create projected PSF raw output");
    if (fwrite(output, sizeof(float), nout, f) != nout) die("short PSF raw output write");
    if (fclose(f)) die("cannot close projected PSF raw output");

    printf("PSF polygon %d/%d: %dx%d -> %dx%d, overlaps=%llu, flux rel=%.3e\n",
           index + 1, njobs, iw, ih, ow, oh,
           (unsigned long long)overlaps, job->flux_rel_error);
    fflush(stdout);
    free(output);
    free(vertices);
    free(input);
}

typedef struct {
    int next;
    pthread_mutex_t mutex;
} Queue;

static void *worker_main(void *opaque) {
    Queue *q = (Queue *)opaque;
    for (;;) {
        pthread_mutex_lock(&q->mutex);
        int i = q->next++;
        pthread_mutex_unlock(&q->mutex);
        if (i >= njobs) break;
        project_one(&jobs[i], i);
    }
    return NULL;
}

static void parse_manifest(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) die("cannot open PSF projection manifest");
    int cap = 16;
    jobs = xcalloc((size_t)cap, sizeof(*jobs));
    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
        if (njobs == cap) {
            cap *= 2;
            Job *tmp = realloc(jobs, (size_t)cap * sizeof(*jobs));
            if (!tmp) die("out of memory expanding jobs");
            jobs = tmp;
            memset(jobs + njobs, 0, (size_t)(cap - njobs) * sizeof(*jobs));
        }
        char *save = NULL;
        char *field[6];
        for (int j = 0; j < 6; ++j) field[j] = strtok_r(j ? NULL : line, "\t\r\n", &save);
        for (int j = 0; j < 6; ++j) if (!field[j]) die("PSF manifest needs input,width,height,vertices,output,label");
        Job *job = &jobs[njobs++];
        if (strlen(field[0]) >= sizeof(job->input) || strlen(field[3]) >= sizeof(job->vertices) ||
            strlen(field[4]) >= sizeof(job->output)) die("PSF manifest path too long");
        strcpy(job->input, field[0]);
        job->width = (int)strtol(field[1], NULL, 10);
        job->height = (int)strtol(field[2], NULL, 10);
        strcpy(job->vertices, field[3]);
        strcpy(job->output, field[4]);
        if (job->width < 1 || job->height < 1) die("invalid PSF input dimensions");
    }
    fclose(f);
    if (njobs < 1) die("empty PSF projection manifest");
}

static void write_results(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) die("cannot create PSF projection results");
    fprintf(f, "exposure,width,height,input_sum,output_sum,flux_rel_error,overlaps\n");
    for (int i = 0; i < njobs; ++i) {
        const Job *j = &jobs[i];
        fprintf(f, "%d,%d,%d,%.17g,%.17g,%.17g,%llu\n",
                i + 1, j->out_width, j->out_height, j->input_sum, j->output_sum,
                j->flux_rel_error, (unsigned long long)j->overlap_count);
    }
    if (fclose(f)) die("cannot close PSF projection results");
}

static int self_test(void) {
    const UpdcPoint quad[4] = {{-0.5,-0.5},{-0.5,0.5},{0.5,0.5},{0.5,-0.5}};
    if (fabsf(updc_polygon_area(quad, 4) - 1.0f) > 2.0e-6f) return 1;
    if (fabsf(updc_overlap_area(quad, 0, 0) - 1.0f) > 2.0e-6f) return 1;
    printf("UPDC PSF projector self-test: OK\n");
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 2 && strcmp(argv[1], "--self-test") == 0) return self_test();
    if (argc != 4) die("usage: updc_psf_projector manifest workers results.csv");
    parse_manifest(argv[1]);
    int workers = atoi(argv[2]);
    if (workers < 1) workers = 1;
    if (workers > njobs) workers = njobs;
    if (workers == 1) {
        for (int i = 0; i < njobs; ++i) project_one(&jobs[i], i);
    } else {
        pthread_t *threads = xmalloc((size_t)workers * sizeof(*threads));
        Queue q = {0};
        if (pthread_mutex_init(&q.mutex, NULL)) die("cannot initialize PSF worker mutex");
        for (int i = 0; i < workers; ++i)
            if (pthread_create(&threads[i], NULL, worker_main, &q)) die("cannot create PSF worker");
        for (int i = 0; i < workers; ++i)
            if (pthread_join(threads[i], NULL)) die("cannot join PSF worker");
        pthread_mutex_destroy(&q.mutex);
        free(threads);
    }
    write_results(argv[3]);
    free(jobs);
    return 0;
}
