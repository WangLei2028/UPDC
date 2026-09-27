/* Signed-x, background-aware multiplicative UPDC with persistent mu freezing.
 *
 * H_k = S_k P_k, where S is exact detector-pixel area integration:
 *   S[p,i] = overlap_area[p,i] / full_mapped_detector_pixel_area[p].
 *
 * Input Level-2 SCI is NOT background-subtracted by the driver.  For exposure k:
 *   y_eff = SCI + poisson_pedestal
 *   raw_mu = H_k x + background_k + poisson_pedestal
 *   ratio  = y_eff/raw_mu while raw_mu > 0; permanently 1 after raw_mu <= 0
 *   x <- x * sum H_k^T[W_k * ratio] / sum H_k^T W_k
 *
 * The pedestal is a constant mean only; it never adds a random noise realization.
 */
#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#define _DARWIN_C_SOURCE 1
#endif
#include <fftw3.h>
#include <math.h>
#include <stdint.h>
#include <stdatomic.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/utsname.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#ifdef UPDC_PARALLEL
#include <pthread.h>
#endif
#include "overlap_geometry.h"

#define EPS_MU 1.0e-12f

/* Versioned, telescope-agnostic LinkMap.  Each exposure independently uses
 * DELTA16 only when every actual mapped detector-pixel stencil fits in a
 * 16-bit target-index delta; otherwise it falls back to ABS32.  No regular
 * stencil assumption participates in the numerical path. */
#define LINKMAP_MAGIC 0x324d4c55u /* "ULM2" on little-endian hosts */
#define LINKMAP_VERSION 2u
#define LINKMAP_ENCODING_ABS32 1u
#define LINKMAP_ENCODING_DELTA16 2u
typedef struct {
    uint32_t magic,version,encoding,header_bytes;
    uint64_t nin,nlinks;
    uint32_t target_width,target_height;
    uint32_t reserved[6];
} LinkMapHeader;
_Static_assert(sizeof(LinkMapHeader)==64,"LinkMapHeader must remain 64 bytes");
typedef struct { uint32_t target; float frac; } TempLink;
typedef struct {
    char science[1024], weight[1024], psf[1024], vertices[1024], links[1100];
    float background, pedestal, mean_offset;
    int psf_width, psf_height;
    FILE *spectrum;
    const fftwf_complex *spectrum_map;
    size_t spectrum_bytes;
} Exposure;
typedef struct {
    int fd; size_t bytes; void *base; uint32_t encoding; uint64_t nlinks;
    const uint32_t *offset;
    const uint32_t *base_target;
    const uint16_t *delta16;
    const uint32_t *target32;
    const float *frac;
} LinkMap;

static size_t align_up_size(size_t value,size_t alignment) {
    return (value + alignment - 1u) & ~(alignment - 1u);
}
static inline uint32_t link_target_at(const LinkMap *m,size_t p,uint32_t z) {
    return m->encoding==LINKMAP_ENCODING_DELTA16 ? m->base_target[p]+(uint32_t)m->delta16[z] : m->target32[z];
}
static inline float link_frac_at(const LinkMap *m,uint32_t z) { return m->frac[z]; }

typedef struct {
    float *blur;
    float *scatter;
    float *fft_image;
    fftwf_complex *fft_freq;
    fftwf_complex *partial_freq;
    double time_forward_ifft,time_detector,time_scatter_fft;
    int owns_buffers;
} ExposureWorkspace;

static Exposure *expv;
static LinkMap *linkmap;
static int W,H,FW,FH,iw,ih,nexp,psf_transposed;
static int fft_workers=1,exposure_workers=1,preprocess_workers=1,fft_plan_measure=1,fft_autotune_enabled=1,fft_retune=0,save_every=0;
static size_t nout,nf,nc,nin;
static float *estimate,*correction,*denominator,*scatter,*coverage,*blur,*science_all,*weight_all,*fft_image;
static ExposureWorkspace *exp_workspace;
static int n_exp_workspace=1;
static uint8_t *mu_frozen;
static float weight_scale;
static float max_sensitivity;
static fftwf_complex *fft_freq,*estimate_freq;
static fftwf_plan forward_plan,inverse_plan;
static atomic_uint_fast64_t fft_forward_count,fft_inverse_count;
static const char *fftw_wisdom_path=NULL;
static char cache_root[4096];
static int persistent_cache_enabled=0;

static void fail(const char *s) { fprintf(stderr,"UPDC_ALL: %s\n",s); exit(EXIT_FAILURE); }
static double wall_seconds(void) {
    struct timespec ts;
    if(clock_gettime(CLOCK_MONOTONIC,&ts))fail("cannot read monotonic clock");
    return (double)ts.tv_sec+1.0e-9*(double)ts.tv_nsec;
}
static void *alloc(size_t n) { void *p=malloc(n); if(!p)fail("out of memory"); return p; }
static void *zalloc(size_t n) { void *p=calloc(1,n); if(!p)fail("out of memory"); return p; }
static int next_fast_size(int n) {
    for(;;++n) {
        int q=n;
        while(q%2==0)q/=2;
        while(q%3==0)q/=3;
        while(q%5==0)q/=5;
        while(q%7==0)q/=7;
        if(q==1)return n;
    }
}
static int cmp_float(const void *a,const void *b) {
    float x=*(const float *)a,y=*(const float *)b; return (x>y)-(x<y);
}
static void set_cache_root_from_manifest(const char *manifest) {
    if(strlen(manifest)>=sizeof(cache_root))fail("manifest path too long for cache root");
    strcpy(cache_root,manifest);
    char *slash=strrchr(cache_root,'/');
#ifdef _WIN32
    char *bslash=strrchr(cache_root,'\\'); if(!slash || (bslash && bslash>slash))slash=bslash;
#endif
    if(slash)*slash='\0'; else strcpy(cache_root,".");
}
static int exact_file_size(const char *path,size_t bytes) {
    struct stat st; return stat(path,&st)==0 && (size_t)st.st_size==bytes;
}
static void cache_path(char *out,size_t cap,const char *name) {
    if(snprintf(out,cap,"%s/%s",cache_root,name)>=(int)cap)fail("cache path too long");
}
static int load_f32_cache(const char *path,float *dest,size_t count) {
    if(!exact_file_size(path,count*sizeof(float)))return 0;
    FILE *f=fopen(path,"rb"); if(!f)return 0;
    const int ok=fread(dest,sizeof(float),count,f)==count && fclose(f)==0;
    return ok;
}
static void save_f32_cache(const char *path,const float *src,size_t count) {
    char tmp[4096]; if(snprintf(tmp,sizeof(tmp),"%s.tmp.%ld",path,(long)getpid())>=(int)sizeof(tmp))fail("cache temp path too long");
    FILE *f=fopen(tmp,"wb"); if(!f)fail("cannot create cache file");
    if(fwrite(src,sizeof(float),count,f)!=count)fail("cache write failed");
    if(fflush(f)||fsync(fileno(f))||fclose(f))fail("cache close failed");
    if(rename(tmp,path))fail("cache rename failed");
}

static void configure_fft_threads(void) {
#ifdef UPDC_FFTW_THREADS
    if(fft_workers>1) {
        if(!fftwf_init_threads())fail("FFTW thread initialization failed");
        fftwf_plan_with_nthreads(fft_workers);
    }
#else
    if(fft_workers>1)fail("binary was built without FFTW thread support");
#endif
}

static int fill_fast_candidates(int minimum,int *out,int cap) {
    int n=0,value=minimum;
    while(n<cap) {
        value=next_fast_size(value);
        out[n++]=value;
        if(value==INT_MAX)break;
        value++;
    }
    return n;
}

static void current_host_signature(char *out,size_t cap) {
    struct utsname u;
    if(uname(&u)==0) {
        if(snprintf(out,cap,"%s|%s|%s|ptr%zu|%s",u.sysname,u.release,u.machine,sizeof(void*)*8,fftwf_version)>=(int)cap)
            fail("host/FFTW signature too long");
    } else {
        if(snprintf(out,cap,"unknown-host|ptr%zu|%s",sizeof(void*)*8,fftwf_version)>=(int)cap)
            fail("host/FFTW signature too long");
    }
}

static int load_fft_autotune_cache(int min_width,int min_height,int *width,int *height,double *seconds) {
    if(!persistent_cache_enabled || fft_retune)return 0;
    char path[4096],name[128];
    if(snprintf(name,sizeof(name),"fft_autotune_v1_w%d_%dx%d.txt",fft_workers,min_width,min_height)>=(int)sizeof(name))fail("FFT autotune cache name too long");
    cache_path(path,sizeof(path),name);
    FILE *f=fopen(path,"r"); if(!f)return 0;
    char host[768]={0},cached_host[768]={0};
    int cw=0,ch=0,cminw=0,cminh=0,cworkers=0;
    double csec=0.0;
    current_host_signature(host,sizeof(host));
    int ok=fscanf(f,"signature=%767[^\n]\nworkers=%d\nmin_width=%d\nmin_height=%d\nwidth=%d\nheight=%d\nroundtrip=%lf\n",
                  cached_host,&cworkers,&cminw,&cminh,&cw,&ch,&csec)==7;
    fclose(f);
    if(!ok || strcmp(host,cached_host)!=0 || cworkers!=fft_workers || cminw!=min_width || cminh!=min_height)return 0;
    if(cw<min_width || ch<min_height)return 0;
    *width=cw;*height=ch;if(seconds)*seconds=csec;
    return 1;
}

static void save_fft_autotune_cache(int min_width,int min_height,int width,int height,double seconds) {
    if(!persistent_cache_enabled)return;
    char path[4096],tmp[4096],host[768],name[128];
    if(snprintf(name,sizeof(name),"fft_autotune_v1_w%d_%dx%d.txt",fft_workers,min_width,min_height)>=(int)sizeof(name))fail("FFT autotune cache name too long");
    cache_path(path,sizeof(path),name);
    current_host_signature(host,sizeof(host));
    if(snprintf(tmp,sizeof(tmp),"%s.tmp.%ld",path,(long)getpid())>=(int)sizeof(tmp))fail("FFT autotune cache temp path too long");
    FILE *f=fopen(tmp,"w"); if(!f)return;
    fprintf(f,"signature=%s\nworkers=%d\nmin_width=%d\nmin_height=%d\nwidth=%d\nheight=%d\nroundtrip=%.9g\n",
            host,fft_workers,min_width,min_height,width,height,seconds);
    if(fflush(f)||fsync(fileno(f))||fclose(f)) {unlink(tmp);return;}
    if(rename(tmp,path))unlink(tmp);
}

static double benchmark_fft_shape(int width,int height,float *image,fftwf_complex *freq) {
    fftwf_plan fwd=fftwf_plan_dft_r2c_2d(height,width,image,freq,FFTW_ESTIMATE);
    fftwf_plan inv=fftwf_plan_dft_c2r_2d(height,width,freq,image,FFTW_ESTIMATE);
    if(!fwd||!inv)fail("FFT autotune planning failed");
    size_t count=(size_t)width*height;
    memset(image,0,count*sizeof(float));
    image[0]=1.0f;
    if(count>1)image[count/2]=0.25f;
    double t0=wall_seconds();
    fftwf_execute(fwd);
    fftwf_execute(inv);
    double elapsed=wall_seconds()-t0;
    fftwf_destroy_plan(fwd);
    fftwf_destroy_plan(inv);
    return elapsed;
}

static void select_fft_shape(int min_width,int min_height) {
    int wc[5],hc[5];
    int nw=fill_fast_candidates(min_width,wc,5);
    int nh=fill_fast_candidates(min_height,hc,5);
    if(!fft_autotune_enabled||nw<1||nh<1) {
        FW=next_fast_size(min_width);
        FH=next_fast_size(min_height);
        return;
    }
    double cached_seconds=0.0;
    if(load_fft_autotune_cache(min_width,min_height,&FW,&FH,&cached_seconds)) {
        printf("FFT autotune persistent-cache hit: %dx%d (roundtrip %.6f s, workers=%d)\n",
               FW,FH,cached_seconds,fft_workers);fflush(stdout);
        return;
    }
    int maxw=wc[nw-1],maxh=hc[nh-1];
    size_t max_nf=(size_t)maxw*maxh;
    size_t max_nc=(size_t)maxh*(maxw/2+1);
    float *bench_image=fftwf_alloc_real(max_nf);
    fftwf_complex *bench_freq=fftwf_alloc_complex(max_nc);
    if(!bench_image||!bench_freq)fail("FFT autotune allocation failed");

    int shapes[16][2],ns=0;
#define ADD_SHAPE(A,B) do { \
        int _w=(A),_h=(B),_dup=0; \
        for(int _i=0;_i<ns;++_i)if(shapes[_i][0]==_w&&shapes[_i][1]==_h)_dup=1; \
        if(!_dup&&ns<(int)(sizeof(shapes)/sizeof(shapes[0]))) {shapes[ns][0]=_w;shapes[ns][1]=_h;++ns;} \
    } while(0)
    int ndiag=nw<nh?nw:nh;
    for(int i=0;i<ndiag;++i)ADD_SHAPE(wc[i],hc[i]);
    if(nw>2)ADD_SHAPE(wc[2],hc[0]);
    if(nh>2)ADD_SHAPE(wc[0],hc[2]);
    if(nw>4)ADD_SHAPE(wc[4],hc[0]);
    if(nh>4)ADD_SHAPE(wc[0],hc[4]);
#undef ADD_SHAPE

    double best=INFINITY; int bestw=wc[0],besth=hc[0];
    printf("FFT autotune: %d candidate shapes, workers=%d\n",ns,fft_workers);fflush(stdout);
    for(int i=0;i<ns;++i) {
        int cw=shapes[i][0],ch=shapes[i][1];
        double elapsed=benchmark_fft_shape(cw,ch,bench_image,bench_freq);
        printf("  FFT candidate %dx%d roundtrip %.6f s\n",cw,ch,elapsed);fflush(stdout);
        if(elapsed<best) {best=elapsed;bestw=cw;besth=ch;}
    }
    FW=bestw;FH=besth;
    printf("FFT autotune selected %dx%d (roundtrip %.6f s)\n",FW,FH,best);fflush(stdout);
    save_fft_autotune_cache(min_width,min_height,FW,FH,best);
    fftwf_free(bench_image);
    fftwf_free(bench_freq);
}

static int host_is_little_endian(void) {
    const uint16_t one=1;
    return *((const uint8_t *)&one)==1;
}

static void read_f32_file_exact(const char *path,float *out,size_t count,const char *label) {
    FILE *f=fopen(path,"rb");
    if(!f) {
        fprintf(stderr,"%s: %s\n",label,path);
        fail("cannot open float32 input");
    }
    if(fread(out,sizeof(float),count,f)!=count) {
        fclose(f);
        fprintf(stderr,"%s: %s\n",label,path);
        fail("short float32 input");
    }
    if(fgetc(f)!=EOF) {
        fclose(f);
        fprintf(stderr,"%s: %s\n",label,path);
        fail("float32 input has unexpected trailing bytes");
    }
    if(fclose(f))fail("cannot close float32 input");
}

#ifdef UPDC_PARALLEL
typedef struct { int worker_id,nthreads; } PreExposureArg;
static void *load_input_worker(void *opaque) {
    PreExposureArg *a=(PreExposureArg *)opaque;
    for(int k=a->worker_id;k<nexp;k+=a->nthreads) {
        read_f32_file_exact(expv[k].science,science_all+(size_t)k*nin,nin,"SCI raw");
        read_f32_file_exact(expv[k].weight,weight_all+(size_t)k*nin,nin,"WHT raw");
    }
    return NULL;
}
#endif
static void load_input_planes(void) {
    science_all=alloc((size_t)nexp*nin*sizeof(float));
    weight_all=alloc((size_t)nexp*nin*sizeof(float));
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        PreExposureArg *args=alloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            args[w]=(PreExposureArg){w,nthreads};
            if(pthread_create(&threads[w],NULL,load_input_worker,&args[w]))fail("cannot create input loader");
        }
        for(int w=0;w<nthreads;++w)if(pthread_join(threads[w],NULL))fail("cannot join input loader");
        free(args);free(threads);
    } else
#endif
    {
        for(int k=0;k<nexp;++k) {
            read_f32_file_exact(expv[k].science,science_all+(size_t)k*nin,nin,"SCI raw");
            read_f32_file_exact(expv[k].weight,weight_all+(size_t)k*nin,nin,"WHT raw");
        }
    }
    printf("SCI/WHT raw planes cached in RAM: %.3f GiB (%d preprocess worker%s)\n",
           (double)(2*(size_t)nexp*nin*sizeof(float))/1073741824.0,
           nthreads,nthreads==1?"":"s");
    fflush(stdout);
}

#ifdef UPDC_PARALLEL
typedef struct { int worker_id,nthreads; size_t active,masked; } WeightNormArg;
static void *weight_norm_worker(void *opaque) {
    WeightNormArg *a=(WeightNormArg *)opaque;
    size_t active=0,masked=0;
    for(int k=a->worker_id;k<nexp;k+=a->nthreads) {
        float *input=science_all+(size_t)k*nin;
        float *weight=weight_all+(size_t)k*nin;
        for(size_t p=0;p<nin;++p) {
            if(!isfinite(input[p])||!isfinite(weight[p])||!(weight[p]>0)) {
                input[p]=0.0f; weight[p]=0.0f; ++masked; continue;
            }
            weight[p]/=weight_scale;
            if(!isfinite(weight[p])||!(weight[p]>0))fail("invalid normalized WHT");
            const float y=input[p]+expv[k].pedestal;
            if(!isfinite(y)||y<0.0f)fail("negative/non-finite y_eff reached solver");
            ++active;
        }
    }
    a->active=active;a->masked=masked;return NULL;
}
#endif
static void determine_weight_scale(void) {
    float *sample=alloc((nexp*(nin/64+1))*sizeof(float)); size_t n=0;
    for(int k=0;k<nexp;++k) {
        const float *weight=weight_all+(size_t)k*nin;
        for(size_t p=0;p<nin;p+=64)
            if(isfinite(weight[p])&&weight[p]>0)sample[n++]=weight[p];
    }
    if(!n)fail("all input WHT pixels are zero");
    qsort(sample,n,sizeof(float),cmp_float); weight_scale=sample[n/2]; free(sample);
    if(!(weight_scale>0)||!isfinite(weight_scale))fail("invalid WHT normalization");

    size_t active=0,masked=0;
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        WeightNormArg *args=zalloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            args[w].worker_id=w;args[w].nthreads=nthreads;
            if(pthread_create(&threads[w],NULL,weight_norm_worker,&args[w]))fail("cannot create WHT normalization worker");
        }
        for(int w=0;w<nthreads;++w) {
            if(pthread_join(threads[w],NULL))fail("cannot join WHT normalization worker");
            active+=args[w].active;masked+=args[w].masked;
        }
        free(args);free(threads);
    } else
#endif
    {
        for(int k=0;k<nexp;++k) {
            float *input=science_all+(size_t)k*nin;
            float *weight=weight_all+(size_t)k*nin;
            for(size_t p=0;p<nin;++p) {
                if(!isfinite(input[p])||!isfinite(weight[p])||!(weight[p]>0)) {
                    input[p]=0.0f;weight[p]=0.0f;++masked;continue;
                }
                weight[p]/=weight_scale;
                if(!isfinite(weight[p])||!(weight[p]>0))fail("invalid normalized WHT");
                const float y=input[p]+expv[k].pedestal;
                if(!isfinite(y)||y<0.0f)fail("negative/non-finite y_eff reached solver");
                ++active;
            }
        }
    }
    printf("WHT median normalization %.9g; normalized once in RAM (%zu active, %zu masked; %d preprocess worker%s)\n",
           weight_scale,active,masked,nthreads,nthreads==1?"":"s"); fflush(stdout);
}

static void parse_manifest(const char *path) {
    FILE *f=fopen(path,"r"); if(!f)fail("cannot open manifest");
    expv=zalloc((size_t)nexp*sizeof(*expv));
    linkmap=zalloc((size_t)nexp*sizeof(*linkmap));
    for(int k=0;k<nexp;++k) {
        char line[8192];
        if(!fgets(line,sizeof(line),f))fail("manifest has fewer rows than expected");
        char *save=NULL,*field[8];
        for(int j=0;j<8;++j)field[j]=strtok_r(j?NULL:line,"\t\r\n",&save);
        for(int j=0;j<8;++j)if(!field[j])
            fail("manifest requires science, weight, psf, psf_width, psf_height, vertices, background, pedestal");
        if(strlen(field[0])>=sizeof(expv[k].science) ||
           strlen(field[1])>=sizeof(expv[k].weight) ||
           strlen(field[2])>=sizeof(expv[k].psf) ||
           strlen(field[5])>=sizeof(expv[k].vertices))fail("manifest path too long");
        strcpy(expv[k].science,field[0]);
        strcpy(expv[k].weight,field[1]);
        strcpy(expv[k].psf,field[2]);
        expv[k].psf_width=(int)strtol(field[3],NULL,10);
        expv[k].psf_height=(int)strtol(field[4],NULL,10);
        strcpy(expv[k].vertices,field[5]);
        expv[k].background=strtof(field[6],NULL);
        expv[k].pedestal=strtof(field[7],NULL);
        expv[k].mean_offset=expv[k].background+expv[k].pedestal;
        if(expv[k].psf_width<1||expv[k].psf_height<1)
            fail("invalid PSF dimensions in manifest");
        if(snprintf(expv[k].links,sizeof(expv[k].links),"%s.links",expv[k].vertices)>=(int)sizeof(expv[k].links))
            fail("link path is too long");
        if(!isfinite(expv[k].background)||!isfinite(expv[k].pedestal)||expv[k].pedestal<0.0f||
           !isfinite(expv[k].mean_offset))
            fail("invalid background/pedestal in manifest");
    }
    char extra[8192];
    if(fgets(extra,sizeof(extra),f))fail("manifest has more rows than expected");
    fclose(f);
}

/* LinkMap v2 layout:
 *   header[64], uint32 offset[nin+1], then either
 *   DELTA16: uint32 base_target[nin], uint16 delta[nlinks], pad-to-4, float frac[nlinks]
 *   ABS32:   uint32 target[nlinks], float frac[nlinks]
 * Offsets remain uint32 for compact hot-loop access; a generic dataset that
 * exceeds UINT32_MAX overlaps is rejected explicitly instead of truncating. */
static size_t linkmap_expected_bytes(const LinkMapHeader *h) {
    if(h->header_bytes!=sizeof(LinkMapHeader))return 0;
    if(h->nin!=nin || h->nlinks>UINT32_MAX)return 0;
    size_t pos=sizeof(LinkMapHeader)+(nin+1)*sizeof(uint32_t);
    if(h->encoding==LINKMAP_ENCODING_DELTA16) {
        pos+=nin*sizeof(uint32_t);
        pos+=(size_t)h->nlinks*sizeof(uint16_t);
        pos=align_up_size(pos,4);
        pos+=(size_t)h->nlinks*sizeof(float);
    } else if(h->encoding==LINKMAP_ENCODING_ABS32) {
        pos+=(size_t)h->nlinks*sizeof(uint32_t);
        pos+=(size_t)h->nlinks*sizeof(float);
    } else return 0;
    return pos;
}

static int read_linkmap_header(const char *path,LinkMapHeader *h,size_t *file_bytes) {
    struct stat st;
    if(stat(path,&st) || st.st_size<(off_t)sizeof(*h))return 0;
    FILE *f=fopen(path,"rb");if(!f)return 0;
    int ok=fread(h,sizeof(*h),1,f)==1;fclose(f);
    if(!ok || h->magic!=LINKMAP_MAGIC || h->version!=LINKMAP_VERSION ||
       h->target_width!=(uint32_t)W || h->target_height!=(uint32_t)H)return 0;
    size_t expected=linkmap_expected_bytes(h);
    if(!expected || (size_t)st.st_size!=expected)return 0;
    if(file_bytes)*file_bytes=expected;
    return 1;
}

static int link_cache_valid(const char *path) {
    if(!persistent_cache_enabled)return 0;
    LinkMapHeader h;return read_linkmap_header(path,&h,NULL);
}

static int is_regular_3x3(const uint32_t t[9]) {
    const uint32_t base=t[0];
    const uint32_t x=base%(uint32_t)W;
    if(x+2u>=(uint32_t)W)return 0;
    for(uint32_t dy=0;dy<3;++dy)for(uint32_t dx=0;dx<3;++dx)
        if(t[dy*3+dx]!=base+dy*(uint32_t)W+dx)return 0;
    return 1;
}

static int mapped_pixel_quad(const double *vertices,int x,int y,UpdcPoint quad[4],
                             float *full_area,int *x0,int *x1,int *y0,int *y1) {
    size_t a=2*((size_t)y*(iw+1)+x),b=a+2*(iw+1);
    quad[0]=(UpdcPoint){vertices[a],vertices[a+1]};quad[1]=(UpdcPoint){vertices[b],vertices[b+1]};
    quad[2]=(UpdcPoint){vertices[b+2],vertices[b+3]};quad[3]=(UpdcPoint){vertices[a+2],vertices[a+3]};
    *full_area=updc_polygon_area(quad,4);
    if(!isfinite(*full_area)||*full_area<=1.0e-12f)return 0;
    double minx=quad[0].x,maxx=minx,miny=quad[0].y,maxy=miny;
    for(int j=1;j<4;++j){if(quad[j].x<minx)minx=quad[j].x;if(quad[j].x>maxx)maxx=quad[j].x;if(quad[j].y<miny)miny=quad[j].y;if(quad[j].y>maxy)maxy=quad[j].y;}
    if(!isfinite(minx)||!isfinite(maxx)||!isfinite(miny)||!isfinite(maxy)||maxx<-.5||maxy<-.5||minx>W-.5||miny>H-.5)return 0;
    *x0=(int)floor(minx+.5);*x1=(int)floor(maxx+.5);*y0=(int)floor(miny+.5);*y1=(int)floor(maxy+.5);
    if(*x0<0)*x0=0;
    if(*y0<0)*y0=0;
    if(*x1>=W)*x1=W-1;
    if(*y1>=H)*y1=H-1;
    return *x0<=*x1 && *y0<=*y1;
}

static void AreaLinkList_O2TAD_stream(int k) {
    if(link_cache_valid(expv[k].links)) {
        LinkMapHeader h;read_linkmap_header(expv[k].links,&h,NULL);
        printf("AreaLink exposure %d: persistent cache hit (%s)\n",k+1,
               h.encoding==LINKMAP_ENCODING_DELTA16?"DELTA16":"ABS32");fflush(stdout);return;
    }
    size_t nv=2*(size_t)(iw+1)*(ih+1);
    double *vertices=alloc(nv*sizeof(double));
    FILE *vf=fopen(expv[k].vertices,"rb");if(!vf)fail("cannot open vertex map");
    if(fread(vertices,sizeof(double),nv,vf)!=nv)fail("short vertex map");
    fclose(vf);
    uint32_t *offset=alloc((nin+1)*sizeof(uint32_t));
    uint32_t *base_target=zalloc(nin*sizeof(uint32_t));
    uint64_t count=0,hist[17]={0},active_pixels=0,span16_pixels=0,nine_pixels=0,regular3_pixels=0;
    uint32_t max_links=0,max_span=0;int all_delta16=1;UpdcPoint quad[4];

    /* Pass 1: exact counts/offsets and encoding capability.  This intentionally
     * repeats geometry in pass 2 to avoid writing a multi-GiB temporary ABS32
     * LinkMap before compression. */
    for(int y=0;y<ih;++y)for(int x=0;x<iw;++x) {
        size_t p=(size_t)y*iw+x;if(count>UINT32_MAX)fail("LinkMap exceeds uint32 offset capacity");offset[p]=(uint32_t)count;
        float full_area;int x0,x1,y0,y1;
        if(!mapped_pixel_quad(vertices,x,y,quad,&full_area,&x0,&x1,&y0,&y1)){hist[0]++;continue;}
        const uint64_t before=count;uint32_t first=0,last=0,nine[9]={0};
        for(int ty=y0;ty<=y1;++ty)for(int tx=x0;tx<=x1;++tx) {
            const float area=updc_overlap_area(quad,tx,ty);if(area<=1.0e-5f)continue;
            if(count>=UINT32_MAX)fail("LinkMap exceeds uint32 offset capacity");
            const uint32_t target=(uint32_t)((size_t)ty*W+tx),local=(uint32_t)(count-before);
            if(count==before)first=target;
            last=target;
            if(local<9)nine[local]=target;
            ++count;
        }
        const uint32_t nlink=(uint32_t)(count-before);hist[nlink<16?nlink:16]++;
        if(nlink){active_pixels++;base_target[p]=first;if(nlink>max_links)max_links=nlink;const uint32_t span=last-first;if(span>max_span)max_span=span;
            if(span<=UINT16_MAX)span16_pixels++;else all_delta16=0;
            if(nlink==9){nine_pixels++;if(is_regular_3x3(nine))regular3_pixels++;}
        }
    }
    offset[nin]=(uint32_t)count;
    const uint32_t encoding=all_delta16?LINKMAP_ENCODING_DELTA16:LINKMAP_ENCODING_ABS32;
    LinkMapHeader h={0};h.magic=LINKMAP_MAGIC;h.version=LINKMAP_VERSION;h.encoding=encoding;h.header_bytes=sizeof(h);
    h.nin=nin;h.nlinks=count;h.target_width=(uint32_t)W;h.target_height=(uint32_t)H;
    const size_t bytes=linkmap_expected_bytes(&h);if(!bytes)fail("cannot size LinkMap");
    char outpath[1400];if(snprintf(outpath,sizeof(outpath),"%s.build.%ld",expv[k].links,(long)getpid())>=(int)sizeof(outpath))fail("LinkMap temp path too long");
    int fd=open(outpath,O_RDWR|O_CREAT|O_TRUNC,0600);if(fd<0)fail("cannot create LinkMap");
    if(ftruncate(fd,(off_t)bytes))fail("cannot size LinkMap file");
    void *mapped=mmap(NULL,bytes,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);if(mapped==MAP_FAILED)fail("cannot mmap writable LinkMap");
    char *base=(char *)mapped;memcpy(base,&h,sizeof(h));size_t pos=sizeof(h);
    uint32_t *out_offset=(uint32_t *)(base+pos);memcpy(out_offset,offset,(nin+1)*sizeof(uint32_t));pos+=(nin+1)*sizeof(uint32_t);
    uint32_t *out_base=NULL,*out_target=NULL;uint16_t *out_delta=NULL;
    if(encoding==LINKMAP_ENCODING_DELTA16){out_base=(uint32_t *)(base+pos);memcpy(out_base,base_target,nin*sizeof(uint32_t));pos+=nin*sizeof(uint32_t);out_delta=(uint16_t *)(base+pos);pos+=(size_t)count*sizeof(uint16_t);pos=align_up_size(pos,4);}
    else {out_target=(uint32_t *)(base+pos);pos+=(size_t)count*sizeof(uint32_t);}
    float *out_frac=(float *)(base+pos);

    /* Pass 2: fill the already-selected portable layout directly. */
    for(int y=0;y<ih;++y)for(int x=0;x<iw;++x) {
        const size_t p=(size_t)y*iw+x;float full_area;int x0,x1,y0,y1;uint32_t z=offset[p];
        if(!mapped_pixel_quad(vertices,x,y,quad,&full_area,&x0,&x1,&y0,&y1))continue;
        for(int ty=y0;ty<=y1;++ty)for(int tx=x0;tx<=x1;++tx) {
            const float area=updc_overlap_area(quad,tx,ty);if(area<=1.0e-5f)continue;
            const uint32_t target=(uint32_t)((size_t)ty*W+tx);if(z>=offset[p+1])fail("LinkMap second-pass count overflow");
            if(encoding==LINKMAP_ENCODING_DELTA16){const uint32_t delta=target-base_target[p];if(delta>UINT16_MAX)fail("DELTA16 LinkMap invariant failed");out_delta[z]=(uint16_t)delta;}
            else out_target[z]=target;
            out_frac[z]=area/full_area;++z;
        }
        if(z!=offset[p+1])fail("LinkMap second-pass count mismatch");
    }
    if(msync(mapped,bytes,MS_SYNC))fail("LinkMap msync failed");
    if(fsync(fd))fail("LinkMap fsync failed");
    if(munmap(mapped,bytes)||close(fd))fail("LinkMap close failed");
    if(rename(outpath,expv[k].links))fail("LinkMap rename failed");
    free(offset);free(base_target);free(vertices);
    printf("AreaLink exposure %d: %llu overlaps; encoding=%s; active=%llu, max-links=%u, max-target-span=%u, span<=u16 %.2f%%; 9-link=%llu regular-3x3=%llu (%.2f%% of 9-link); hist 0:%llu 4:%llu 5:%llu 6:%llu 7:%llu 8:%llu 9:%llu >=16:%llu\n",
           k+1,(unsigned long long)count,encoding==LINKMAP_ENCODING_DELTA16?"DELTA16":"ABS32",(unsigned long long)active_pixels,max_links,max_span,
           active_pixels?100.0*(double)span16_pixels/(double)active_pixels:0.0,(unsigned long long)nine_pixels,(unsigned long long)regular3_pixels,
           nine_pixels?100.0*(double)regular3_pixels/(double)nine_pixels:0.0,(unsigned long long)hist[0],(unsigned long long)hist[4],
           (unsigned long long)hist[5],(unsigned long long)hist[6],(unsigned long long)hist[7],(unsigned long long)hist[8],
           (unsigned long long)hist[9],(unsigned long long)hist[16]);fflush(stdout);
}

#ifdef UPDC_PARALLEL
typedef struct { int next; pthread_mutex_t mutex; } GeometryQueue;
static void *geometry_worker(void *opaque) {
    GeometryQueue *q=(GeometryQueue *)opaque;
    for(;;) {
        pthread_mutex_lock(&q->mutex);
        int k=q->next++;
        pthread_mutex_unlock(&q->mutex);
        if(k>=nexp)break;
        AreaLinkList_O2TAD_stream(k);
    }
    return NULL;
}
#endif

static void build_area_links(void) {
#ifdef UPDC_PARALLEL
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        GeometryQueue q={0};
        if(pthread_mutex_init(&q.mutex,NULL))fail("cannot initialize geometry mutex");
        for(int i=0;i<nthreads;++i)
            if(pthread_create(&threads[i],NULL,geometry_worker,&q))fail("cannot create geometry worker");
        for(int i=0;i<nthreads;++i)
            if(pthread_join(threads[i],NULL))fail("cannot join geometry worker");
        pthread_mutex_destroy(&q.mutex);
        free(threads);
        return;
    }
#endif
    for(int k=0;k<nexp;++k)AreaLinkList_O2TAD_stream(k);
}

static LinkMap open_links(int k) {
    LinkMap m={0};LinkMapHeader h;size_t expected=0;
    if(!read_linkmap_header(expv[k].links,&h,&expected))fail("invalid or obsolete LinkMap cache");
    m.fd=open(expv[k].links,O_RDONLY);if(m.fd<0)fail("cannot open area links");m.bytes=expected;m.encoding=h.encoding;m.nlinks=h.nlinks;
    m.base=mmap(NULL,m.bytes,PROT_READ,MAP_PRIVATE,m.fd,0);if(m.base==MAP_FAILED)fail("cannot map area links");
    const char *base=(const char *)m.base;size_t pos=sizeof(LinkMapHeader);
    m.offset=(const uint32_t *)(base+pos);pos+=(nin+1)*sizeof(uint32_t);
    if(m.offset[nin]!=(uint32_t)h.nlinks)fail("LinkMap offset/link-count mismatch");
    if(m.encoding==LINKMAP_ENCODING_DELTA16) {
        m.base_target=(const uint32_t *)(base+pos);pos+=nin*sizeof(uint32_t);
        m.delta16=(const uint16_t *)(base+pos);pos+=(size_t)h.nlinks*sizeof(uint16_t);pos=align_up_size(pos,4);
        m.frac=(const float *)(base+pos);
    } else {
        m.target32=(const uint32_t *)(base+pos);pos+=(size_t)h.nlinks*sizeof(uint32_t);
        m.frac=(const float *)(base+pos);
    }
    return m;
}
static void close_links(LinkMap *m){munmap(m->base,m->bytes);close(m->fd);}

static uint64_t prefault_region(const void *base,size_t bytes) {
    if(!base || !bytes)return 0;
#ifdef MADV_WILLNEED
    (void)madvise((void *)base,bytes,MADV_WILLNEED);
#endif
#ifdef MADV_SEQUENTIAL
    (void)madvise((void *)base,bytes,MADV_SEQUENTIAL);
#endif
    long ps=sysconf(_SC_PAGESIZE);
    size_t page=(ps>0)?(size_t)ps:4096u;
    const volatile unsigned char *p=(const volatile unsigned char *)base;
    uint64_t checksum=0;
    for(size_t off=0;off<bytes;off+=page)checksum+=p[off];
    checksum+=p[bytes-1];
    return checksum;
}

#ifdef UPDC_PARALLEL
typedef struct { int worker_id,nthreads; uint64_t checksum; size_t bytes; } PrewarmArg;
static void *prewarm_worker(void *opaque) {
    PrewarmArg *a=(PrewarmArg *)opaque;
    uint64_t checksum=0;size_t bytes=0;
    for(int k=a->worker_id;k<nexp;k+=a->nthreads) {
        checksum+=prefault_region(linkmap[k].base,linkmap[k].bytes);bytes+=linkmap[k].bytes;
        checksum+=prefault_region(expv[k].spectrum_map,expv[k].spectrum_bytes);bytes+=expv[k].spectrum_bytes;
    }
    a->checksum=checksum;a->bytes=bytes;return NULL;
}
#endif
static double prewarm_mmaps(void) {
    double t0=wall_seconds();uint64_t checksum=0;size_t bytes=0;
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        PrewarmArg *args=zalloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            args[w].worker_id=w;args[w].nthreads=nthreads;
            if(pthread_create(&threads[w],NULL,prewarm_worker,&args[w]))fail("cannot create mmap prewarm worker");
        }
        for(int w=0;w<nthreads;++w) {
            if(pthread_join(threads[w],NULL))fail("cannot join mmap prewarm worker");
            checksum+=args[w].checksum;bytes+=args[w].bytes;
        }
        free(args);free(threads);
    } else
#endif
    {
        for(int k=0;k<nexp;++k) {
            checksum+=prefault_region(linkmap[k].base,linkmap[k].bytes);bytes+=linkmap[k].bytes;
            checksum+=prefault_region(expv[k].spectrum_map,expv[k].spectrum_bytes);bytes+=expv[k].spectrum_bytes;
        }
    }
    double elapsed=wall_seconds()-t0;
    printf("Prewarm mmap working set: %.3f GiB in %.3f s (checksum=%llu; %d preprocess worker%s)\n",
           (double)bytes/1073741824.0,elapsed,(unsigned long long)checksum,
           nthreads,nthreads==1?"":"s");fflush(stdout);
    return elapsed;
}

/* FFT helpers.  The estimate spectrum is computed once per UPDC iteration and
 * reused across every exposure.  PSF spectra are file-backed mmaps, avoiding a
 * ~100 MB fread/copy for every convolution on full mosaics. */
static void pad_target_image(const float *src) {
    memset(fft_image,0,nf*sizeof(float));
    for(int y=0;y<H;++y)
        memcpy(fft_image+(size_t)y*FW,src+(size_t)y*W,(size_t)W*sizeof(float));
}
static void forward_target(const float *src,fftwf_complex *dest) {
    pad_target_image(src);
    fftwf_execute_dft_r2c(forward_plan,fft_image,dest);
    atomic_fetch_add_explicit(&fft_forward_count,1,memory_order_relaxed);
}
static void inverse_crop(float *dest) {
    fftwf_execute_dft_c2r(inverse_plan,fft_freq,fft_image);
    atomic_fetch_add_explicit(&fft_inverse_count,1,memory_order_relaxed);
    const float scale=1.0f/(float)nf;
    for(int y=0;y<H;++y)for(int x=0;x<W;++x)
        dest[(size_t)y*W+x]=fft_image[(size_t)y*FW+x]*scale;
}
static void prepare_estimate_spectrum(void) {
    forward_target(estimate,estimate_freq);
}
static void ws_pad_target(ExposureWorkspace *ws,const float *src) {
    memset(ws->fft_image,0,nf*sizeof(float));
    for(int y=0;y<H;++y)
        memcpy(ws->fft_image+(size_t)y*FW,src+(size_t)y*W,(size_t)W*sizeof(float));
}
static void ws_multiply_spectrum(ExposureWorkspace *ws,const fftwf_complex *src,
                                 const fftwf_complex *kernel,int adjoint) {
    for(size_t p=0;p<nc;++p) {
        float a=src[p][0],b=src[p][1],c=kernel[p][0],d=kernel[p][1];
        ws->fft_freq[p][0]=adjoint?a*c+b*d:a*c-b*d;
        ws->fft_freq[p][1]=adjoint?b*c-a*d:a*d+b*c;
    }
}
static void ws_inverse_crop(ExposureWorkspace *ws,float *dest) {
    fftwf_execute_dft_c2r(inverse_plan,ws->fft_freq,ws->fft_image);
    atomic_fetch_add_explicit(&fft_inverse_count,1,memory_order_relaxed);
    const float scale=1.0f/(float)nf;
    for(int y=0;y<H;++y)for(int x=0;x<W;++x)
        dest[(size_t)y*W+x]=ws->fft_image[(size_t)y*FW+x]*scale;
}
static void ws_forward_target(ExposureWorkspace *ws,const float *src) {
    ws_pad_target(ws,src);
    fftwf_execute_dft_r2c(forward_plan,ws->fft_image,ws->fft_freq);
    atomic_fetch_add_explicit(&fft_forward_count,1,memory_order_relaxed);
}
static void ws_convolve_estimate(ExposureWorkspace *ws,float *dest,const fftwf_complex *kernel) {
    ws_multiply_spectrum(ws,estimate_freq,kernel,0);
    ws_inverse_crop(ws,dest);
}
static float compensated_psf_sum(const float *psf,size_t count) {
    /* Neumaier-style compensated summation, entirely in float32. */
    float sum=0.0f,comp=0.0f;
    for(size_t p=0;p<count;++p) {
        const float v=psf[p];
        if(!isfinite(v)||v<0.0f)fail("invalid PSF pixel");
        const float t=sum+v;
        if(fabsf(sum)>=fabsf(v))comp+=(sum-t)+v;
        else comp+=(v-t)+sum;
        sum=t;
    }
    return sum+comp;
}

static void prepare_psf_with_buffers(int k,float *psf,float *local_fft_image,fftwf_complex *local_fft_freq) {
    char path[4096],name[256];
    if(snprintf(name,sizeof(name),"psf_spectrum_%02d_%dx%d_t%d.c64",k+1,FW,FH,psf_transposed)>=(int)sizeof(name))
        fail("PSF spectrum cache name too long");
    cache_path(path,sizeof(path),name);
    expv[k].spectrum_bytes=nc*sizeof(fftwf_complex);
    if(persistent_cache_enabled && exact_file_size(path,expv[k].spectrum_bytes)) {
        expv[k].spectrum=fopen(path,"rb"); if(!expv[k].spectrum)fail("cannot open cached PSF spectrum");
        const int fd=fileno(expv[k].spectrum);
        void *mapped=mmap(NULL,expv[k].spectrum_bytes,PROT_READ,MAP_PRIVATE,fd,0);
        if(mapped==MAP_FAILED)fail("cannot mmap cached PSF spectrum");
        expv[k].spectrum_map=(const fftwf_complex *)mapped;
        printf("PSF %d size %dx%d spectrum persistent-cache hit (%.1f MiB)\n",
               k+1,expv[k].psf_width,expv[k].psf_height,(double)expv[k].spectrum_bytes/1048576.0);fflush(stdout);
        return;
    }
    const int pw=expv[k].psf_width,ph=expv[k].psf_height;
    const size_t count=(size_t)pw*ph;
    read_f32_file_exact(expv[k].psf,psf,count,"PSF raw");
    const float sum=compensated_psf_sum(psf,count);
    if(!(sum>0.0f)||!isfinite(sum))fail("PSF sum is zero/non-finite");
    const float inv_sum=1.0f/sum;
    memset(local_fft_image,0,nf*sizeof(float));
    const int forward_width=psf_transposed?ph:pw;
    const int forward_height=psf_transposed?pw:ph;
    for(int y=0;y<forward_height;++y)for(int x=0;x<forward_width;++x) {
        const int yy=(y-forward_height/2+FH)%FH;
        const int xx=(x-forward_width/2+FW)%FW;
        const size_t source_index=psf_transposed?(size_t)x*pw+y:(size_t)y*pw+x;
        local_fft_image[(size_t)yy*FW+xx]=psf[source_index]*inv_sum;
    }
    fftwf_execute_dft_r2c(forward_plan,local_fft_image,local_fft_freq);
    atomic_fetch_add_explicit(&fft_forward_count,1,memory_order_relaxed);
    if(persistent_cache_enabled) {
        char tmp[4096]; if(snprintf(tmp,sizeof(tmp),"%s.tmp.%ld",path,(long)getpid())>=(int)sizeof(tmp))fail("PSF cache temp path too long");
        FILE *out=fopen(tmp,"wb"); if(!out)fail("cannot create PSF spectrum cache");
        if(fwrite(local_fft_freq,sizeof(fftwf_complex),nc,out)!=nc)fail("PSF cache write failed");
        if(fflush(out)||fsync(fileno(out))||fclose(out))fail("PSF cache close failed");
        if(rename(tmp,path))fail("PSF cache rename failed");
        expv[k].spectrum=fopen(path,"rb"); if(!expv[k].spectrum)fail("cannot reopen PSF spectrum cache");
    } else {
        expv[k].spectrum=tmpfile(); if(!expv[k].spectrum)fail("cannot create temporary PSF spectrum");
        if(fwrite(local_fft_freq,sizeof(fftwf_complex),nc,expv[k].spectrum)!=nc)fail("temporary PSF spectrum write failed");
        if(fflush(expv[k].spectrum)||fsync(fileno(expv[k].spectrum)))fail("temporary PSF spectrum flush failed");
    }
    const int fd=fileno(expv[k].spectrum);
    void *mapped=mmap(NULL,expv[k].spectrum_bytes,PROT_READ,MAP_PRIVATE,fd,0);
    if(mapped==MAP_FAILED)fail("cannot mmap PSF spectrum");
    expv[k].spectrum_map=(const fftwf_complex *)mapped;
    printf("PSF %d size %dx%d normalized sum %.9g; spectrum %s %.1f MiB\n",
           k+1,pw,ph,(double)sum,persistent_cache_enabled?"persistent-cache write":"temporary mmap",
           (double)expv[k].spectrum_bytes/1048576.0); fflush(stdout);
}

#ifdef UPDC_PARALLEL
typedef struct { int worker_id,nthreads,max_psf_width,max_psf_height; } PsfPrepArg;
static void *psf_prep_worker(void *opaque) {
    PsfPrepArg *a=(PsfPrepArg *)opaque;
    float *psf=alloc((size_t)a->max_psf_width*a->max_psf_height*sizeof(float));
    float *image=fftwf_alloc_real(nf);
    fftwf_complex *freq=fftwf_alloc_complex(nc);
    if(!image||!freq)fail("cannot allocate parallel PSF FFT workspace");
    for(int k=a->worker_id;k<nexp;k+=a->nthreads)
        prepare_psf_with_buffers(k,psf,image,freq);
    free(psf);fftwf_free(image);fftwf_free(freq);
    return NULL;
}
#endif
static void prepare_all_psfs(int max_psf_width,int max_psf_height) {
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        PsfPrepArg *args=alloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            args[w]=(PsfPrepArg){w,nthreads,max_psf_width,max_psf_height};
            if(pthread_create(&threads[w],NULL,psf_prep_worker,&args[w]))fail("cannot create PSF prep worker");
        }
        for(int w=0;w<nthreads;++w)if(pthread_join(threads[w],NULL))fail("cannot join PSF prep worker");
        free(args);free(threads);return;
    }
#endif
    float *psf=alloc((size_t)max_psf_width*max_psf_height*sizeof(float));
    for(int k=0;k<nexp;++k)prepare_psf_with_buffers(k,psf,fft_image,fft_freq);
    free(psf);
}

/* Initial signed estimate: weighted drizzle of raw SCI minus the correspondingly
 * weighted scalar background, without positivity clipping. */
#ifdef UPDC_PARALLEL
typedef struct {
    int worker_id,nthreads;
    float *coverage_part,*data_part,*bkg_part;
} DrizzleArg;
static void *drizzle_worker(void *opaque) {
    DrizzleArg *a=(DrizzleArg *)opaque;
    for(int k=a->worker_id;k<nexp;k+=a->nthreads) {
        const float *input=science_all+(size_t)k*nin;
        const float *weight=weight_all+(size_t)k*nin;
        const LinkMap *m=&linkmap[k];
        for(size_t p=0;p<nin;++p) {
            if(!(weight[p]>0)||m->offset[p]==m->offset[p+1])continue;
            const float wp=weight[p];
            for(uint32_t z=m->offset[p];z<m->offset[p+1];++z) {
                const uint32_t target=link_target_at(m,p,z);const float q=link_frac_at(m,z)*wp;
                a->coverage_part[target]+=q;
                a->data_part[target]=fmaf(q,input[p],a->data_part[target]);
                a->bkg_part[target]=fmaf(q,expv[k].background,a->bkg_part[target]);
            }
        }
    }
    return NULL;
}
#endif
static void SNA_HRA_stream(void) {
    int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
    float signed_sum=0.0f; size_t nnegative=0;
    char estimate_cache[4096],coverage_cache[4096];
    cache_path(estimate_cache,sizeof(estimate_cache),"initial_estimate.f32");
    cache_path(coverage_cache,sizeof(coverage_cache),"initial_coverage.f32");
    if(persistent_cache_enabled && load_f32_cache(estimate_cache,estimate,nout) && load_f32_cache(coverage_cache,coverage,nout)) {
        for(size_t t=0;t<nout;++t) { signed_sum+=estimate[t]; if(estimate[t]<0)nnegative++; }
        printf("Initial drizzle persistent-cache hit: signed sum %.9g; negative pixels %zu/%zu\n",
               (double)signed_sum,nnegative,nout);fflush(stdout);
        return;
    }
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        DrizzleArg *args=zalloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            args[w].worker_id=w;args[w].nthreads=nthreads;
            args[w].coverage_part=zalloc(nout*sizeof(float));
            args[w].data_part=zalloc(nout*sizeof(float));
            args[w].bkg_part=zalloc(nout*sizeof(float));
            if(pthread_create(&threads[w],NULL,drizzle_worker,&args[w]))fail("cannot create drizzle worker");
        }
        for(int w=0;w<nthreads;++w)if(pthread_join(threads[w],NULL))fail("cannot join drizzle worker");
        for(size_t t=0;t<nout;++t) {
            float c=0.0f,d=0.0f,b=0.0f;
            for(int w=0;w<nthreads;++w) {
                c+=args[w].coverage_part[t];d+=args[w].data_part[t];b+=args[w].bkg_part[t];
            }
            coverage[t]=c;
            if(c>1e-12f) {
                estimate[t]=d/c-b/c;
                if(!isfinite(estimate[t]))fail("non-finite signed initial estimate");
                signed_sum+=estimate[t];if(estimate[t]<0)nnegative++;
            } else estimate[t]=0.0f;
        }
        for(int w=0;w<nthreads;++w) {
            free(args[w].coverage_part);free(args[w].data_part);free(args[w].bkg_part);
        }
        free(args);free(threads);
    } else
#endif
    {
        float *data_sum=zalloc(nout*sizeof(float));
        float *bkg_sum=zalloc(nout*sizeof(float));
        for(int k=0;k<nexp;++k) {
            const float *input=science_all+(size_t)k*nin;
            const float *weight=weight_all+(size_t)k*nin;
            const LinkMap *m=&linkmap[k];
            for(size_t p=0;p<nin;++p) {
                if(!(weight[p]>0)||m->offset[p]==m->offset[p+1])continue;
                float wp=weight[p];
                for(uint32_t z=m->offset[p];z<m->offset[p+1];++z) {
                    const uint32_t target=link_target_at(m,p,z);float q=link_frac_at(m,z)*wp;
                    coverage[target]+=q;
                    data_sum[target]=fmaf(q,input[p],data_sum[target]);
                    bkg_sum[target]=fmaf(q,expv[k].background,bkg_sum[target]);
                }
            }
        }
        for(size_t t=0;t<nout;++t) {
            if(coverage[t]>1e-12f) {
                estimate[t]=data_sum[t]/coverage[t]-bkg_sum[t]/coverage[t];
                if(!isfinite(estimate[t]))fail("non-finite signed initial estimate");
                signed_sum+=estimate[t];if(estimate[t]<0)nnegative++;
            } else estimate[t]=0;
        }
        free(data_sum);free(bkg_sum);
    }
    if(persistent_cache_enabled) {
        save_f32_cache(estimate_cache,estimate,nout);
        save_f32_cache(coverage_cache,coverage,nout);
    }
    for(int k=0;k<nexp;++k)
        printf("Initial drizzle exposure %d: background %.9g, pedestal %.9g\n",
               k+1,expv[k].background,expv[k].pedestal);
    printf("Initial signed estimate sum %.9g; negative pixels %zu/%zu (%d preprocess worker%s)\n",
           signed_sum,nnegative,nout,nthreads,nthreads==1?"":"s"); fflush(stdout);
}

/* D = sum_k P_k^T S_k^T W_k.  Linearity lets all adjoint spectra be
 * accumulated first, followed by one inverse FFT instead of one inverse per
 * exposure. */
static void accumulate_adjoint_spectrum(fftwf_complex *dest,const fftwf_complex *src,
                                        const fftwf_complex *kernel) {
    for(size_t p=0;p<nc;++p) {
        const float a=src[p][0],b=src[p][1],c=kernel[p][0],d=kernel[p][1];
        dest[p][0]+=a*c+b*d;
        dest[p][1]+=b*c-a*d;
    }
}

#ifdef UPDC_PARALLEL
typedef struct { int worker_id,nthreads; fftwf_complex *partial_freq; } SensitivityArg;
static void *sensitivity_worker(void *opaque) {
    SensitivityArg *a=(SensitivityArg *)opaque;
    ExposureWorkspace ws={0};
    ws.scatter=zalloc(nout*sizeof(float));
    ws.fft_image=fftwf_alloc_real(nf);
    ws.fft_freq=fftwf_alloc_complex(nc);
    if(!ws.fft_image||!ws.fft_freq)fail("cannot allocate sensitivity FFT workspace");
    memset(a->partial_freq,0,nc*sizeof(fftwf_complex));
    for(int k=a->worker_id;k<nexp;k+=a->nthreads) {
        const float *weight=weight_all+(size_t)k*nin;
        const LinkMap *m=&linkmap[k];
        memset(ws.scatter,0,nout*sizeof(float));
        for(size_t p=0;p<nin;++p) {
            if(!(weight[p]>0)||m->offset[p]==m->offset[p+1])continue;
            const float wp=weight[p];
            for(uint32_t z=m->offset[p];z<m->offset[p+1];++z) {
                const uint32_t target=link_target_at(m,p,z);
                ws.scatter[target]=fmaf(link_frac_at(m,z),wp,ws.scatter[target]);
            }
        }
        ws_forward_target(&ws,ws.scatter);
        accumulate_adjoint_spectrum(a->partial_freq,ws.fft_freq,expv[k].spectrum_map);
    }
    free(ws.scatter);fftwf_free(ws.fft_image);fftwf_free(ws.fft_freq);
    return NULL;
}
#endif
static void precompute_denominator(void) {
    char name[256],path[4096];
    if(snprintf(name,sizeof(name),"denominator_%dx%d.f32",FW,FH)>=(int)sizeof(name))fail("denominator cache name too long");
    cache_path(path,sizeof(path),name);
    if(persistent_cache_enabled && load_f32_cache(path,denominator,nout)) {
        max_sensitivity=0.0f;
        for(size_t t=0;t<nout;++t)if(denominator[t]>max_sensitivity)max_sensitivity=denominator[t];
        if(!(max_sensitivity>0)||!isfinite(max_sensitivity))fail("invalid cached H^T W sensitivity");
        printf("H^T W sensitivity persistent-cache hit (%dx%d FFT)\n",FW,FH);fflush(stdout);
        return;
    }
    const int nthreads=preprocess_workers<nexp?preprocess_workers:nexp;
    fftwf_complex **partials=zalloc((size_t)nthreads*sizeof(*partials));
#ifdef UPDC_PARALLEL
    if(nthreads>1) {
        pthread_t *threads=alloc((size_t)nthreads*sizeof(*threads));
        SensitivityArg *args=zalloc((size_t)nthreads*sizeof(*args));
        for(int w=0;w<nthreads;++w) {
            partials[w]=fftwf_alloc_complex(nc); if(!partials[w])fail("cannot allocate sensitivity partial spectrum");
            args[w].worker_id=w;args[w].nthreads=nthreads;args[w].partial_freq=partials[w];
            if(pthread_create(&threads[w],NULL,sensitivity_worker,&args[w]))fail("cannot create sensitivity worker");
        }
        for(int w=0;w<nthreads;++w)if(pthread_join(threads[w],NULL))fail("cannot join sensitivity worker");
        free(args);free(threads);
    } else
#endif
    {
        partials[0]=fftwf_alloc_complex(nc);if(!partials[0])fail("cannot allocate sensitivity spectrum");
        memset(partials[0],0,nc*sizeof(fftwf_complex));
        ExposureWorkspace ws={0}; ws.scatter=zalloc(nout*sizeof(float)); ws.fft_image=fftwf_alloc_real(nf); ws.fft_freq=fftwf_alloc_complex(nc);
        if(!ws.fft_image||!ws.fft_freq)fail("cannot allocate sensitivity FFT workspace");
        for(int k=0;k<nexp;++k) {
            const float *weight=weight_all+(size_t)k*nin; const LinkMap *m=&linkmap[k];
            memset(ws.scatter,0,nout*sizeof(float));
            for(size_t p=0;p<nin;++p) {
                if(!(weight[p]>0)||m->offset[p]==m->offset[p+1])continue;
                const float wp=weight[p];
                for(uint32_t z=m->offset[p];z<m->offset[p+1];++z) {
                    const uint32_t target=link_target_at(m,p,z);ws.scatter[target]=fmaf(link_frac_at(m,z),wp,ws.scatter[target]);
                }
            }
            ws_forward_target(&ws,ws.scatter);
            accumulate_adjoint_spectrum(partials[0],ws.fft_freq,expv[k].spectrum_map);
        }
        free(ws.scatter);fftwf_free(ws.fft_image);fftwf_free(ws.fft_freq);
    }
    for(size_t p=0;p<nc;++p) {
        float re=0.0f,im=0.0f;
        for(int w=0;w<nthreads;++w){re+=partials[w][p][0];im+=partials[w][p][1];}
        fft_freq[p][0]=re;fft_freq[p][1]=im;
    }
    inverse_crop(denominator); /* exactly one adjoint inverse FFT */
    for(int w=0;w<nthreads;++w)fftwf_free(partials[w]);
    free(partials);
    max_sensitivity=0.0f;
    for(size_t t=0;t<nout;++t)if(denominator[t]>max_sensitivity)max_sensitivity=denominator[t];
    if(!(max_sensitivity>0)||!isfinite(max_sensitivity))fail("invalid H^T W sensitivity");
    if(persistent_cache_enabled)save_f32_cache(path,denominator,nout);
    printf("H^T W sensitivity prepared with %d preprocess worker%s and one combined inverse FFT\n",
           nthreads,nthreads==1?"":"s");fflush(stdout);
}

/* Iteration hot path.  Forward prediction, persistent FreezeMu, and scatter
 * construction are fused while each detector pixel's short LinkMap stencil is
 * still hot in cache.  Each worker accumulates adjoint products in frequency
 * space; the main thread performs one deterministic frequency reduction and one
 * inverse FFT for the full correction. */
static void init_exposure_workspaces(void) {
    n_exp_workspace=exposure_workers<nexp?exposure_workers:nexp;
    if(n_exp_workspace<1)n_exp_workspace=1;
    exp_workspace=zalloc((size_t)n_exp_workspace*sizeof(*exp_workspace));
    exp_workspace[0].blur=blur;
    exp_workspace[0].scatter=scatter;
    exp_workspace[0].fft_image=fft_image;
    exp_workspace[0].fft_freq=fft_freq;
    exp_workspace[0].partial_freq=fftwf_alloc_complex(nc);
    if(!exp_workspace[0].partial_freq)fail("cannot allocate exposure partial spectrum");
    exp_workspace[0].owns_buffers=0;
    size_t extra_bytes=nc*sizeof(fftwf_complex);
    for(int w=1;w<n_exp_workspace;++w) {
        ExposureWorkspace *ws=&exp_workspace[w];
        ws->blur=alloc(nout*sizeof(float));
        ws->scatter=alloc(nout*sizeof(float));
        ws->fft_image=fftwf_alloc_real(nf);
        ws->fft_freq=fftwf_alloc_complex(nc);
        ws->partial_freq=fftwf_alloc_complex(nc);
        if(!ws->fft_image||!ws->fft_freq||!ws->partial_freq)fail("cannot allocate exposure worker FFT workspace");
        memset(ws->blur,0,nout*sizeof(float));memset(ws->scatter,0,nout*sizeof(float));
        memset(ws->fft_image,0,nf*sizeof(float));memset(ws->fft_freq,0,nc*sizeof(fftwf_complex));
        memset(ws->partial_freq,0,nc*sizeof(fftwf_complex));
        ws->owns_buffers=1;
        extra_bytes+=2*nout*sizeof(float)+nf*sizeof(float)+2*nc*sizeof(fftwf_complex);
    }
    printf("Exposure pipeline: %d persistent worker%s x %d FFT thread%s; extra workspace %.3f GiB\n",
           n_exp_workspace,n_exp_workspace==1?"":"s",fft_workers,fft_workers==1?"":"s",
           (double)extra_bytes/1073741824.0);fflush(stdout);
    if(n_exp_workspace>1 && fft_workers>1) {
        printf("WARNING: nested parallelism may use up to %d compute threads during concurrent FFTs; benchmark worker combinations.\n",
               n_exp_workspace*fft_workers);fflush(stdout);
    }
}

static void free_exposure_workspaces(void) {
    if(!exp_workspace)return;
    for(int w=0;w<n_exp_workspace;++w) {
        ExposureWorkspace *ws=&exp_workspace[w];
        if(ws->partial_freq)fftwf_free(ws->partial_freq);
        if(w>0) {
            free(ws->blur);free(ws->scatter);fftwf_free(ws->fft_image);fftwf_free(ws->fft_freq);
        }
    }
    free(exp_workspace);exp_workspace=NULL;
}

static void fused_exposure(int k,ExposureWorkspace *ws,size_t *new_out,size_t *total_out,float *min_out) {
    const float *input=science_all+(size_t)k*nin;
    const float *weight=weight_all+(size_t)k*nin;
    uint8_t *frozen=mu_frozen+(size_t)k*nin;
    const LinkMap *m=&linkmap[k];
    double t0=wall_seconds();
    ws_convolve_estimate(ws,ws->blur,expv[k].spectrum_map);
    ws->time_forward_ifft+=wall_seconds()-t0;

    t0=wall_seconds();
    memset(ws->scatter,0,nout*sizeof(float));
    size_t nnew=0,ntotal=0;float min_raw_mu=INFINITY;
    for(size_t p=0;p<nin;++p) {
        if(!(weight[p]>0)||m->offset[p]==m->offset[p+1])continue;
        float pred=0.0f;
        const uint32_t z0=m->offset[p],z1=m->offset[p+1];
        for(uint32_t z=z0;z<z1;++z) {
            const uint32_t target=link_target_at(m,p,z);pred=fmaf(link_frac_at(m,z),ws->blur[target],pred);
        }
        const float raw_mu=pred+expv[k].mean_offset;
        if(!isfinite(raw_mu))fail("non-finite forward mean");
        if(raw_mu<min_raw_mu)min_raw_mu=raw_mu;
        if(!frozen[p] && !(raw_mu>0.0f)) { frozen[p]=1; ++nnew; }
        float ratio;
        if(frozen[p]) { ++ntotal; ratio=1.0f; }
        else {
            const float y=input[p]+expv[k].pedestal;
            ratio=y/raw_mu;
            if(!isfinite(ratio))fail("non-finite y/raw_mu ratio");
        }
        const float wr=weight[p]*ratio;
        /* Reuse the same tiny contiguous stencil immediately while its links are
         * still resident in L1/L2 cache. */
        for(uint32_t z=z0;z<z1;++z) {
            const uint32_t target=link_target_at(m,p,z);
            ws->scatter[target]=fmaf(link_frac_at(m,z),wr,ws->scatter[target]);
        }
    }
    ws->time_detector+=wall_seconds()-t0;
    *new_out=nnew;*total_out=ntotal;*min_out=min_raw_mu;

    t0=wall_seconds();
    ws_forward_target(ws,ws->scatter);
    accumulate_adjoint_spectrum(ws->partial_freq,ws->fft_freq,expv[k].spectrum_map);
    ws->time_scatter_fft+=wall_seconds()-t0;
}

#ifdef UPDC_PARALLEL
typedef struct IterWorkerArg IterWorkerArg;
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t work_cv,done_cv;
    pthread_t *threads;
    IterWorkerArg *args;
    int generation,stop,pending;
    size_t *new_frozen,*total_frozen;
    float *min_raw_mu;
} IterPool;
struct IterWorkerArg { int worker_id,seen_generation; IterPool *pool; };
static IterPool iter_pool;

static void *persistent_iteration_worker(void *opaque) {
    IterWorkerArg *arg=(IterWorkerArg *)opaque;IterPool *pool=arg->pool;
    for(;;) {
        pthread_mutex_lock(&pool->mutex);
        while(!pool->stop && arg->seen_generation==pool->generation)
            pthread_cond_wait(&pool->work_cv,&pool->mutex);
        if(pool->stop){pthread_mutex_unlock(&pool->mutex);break;}
        arg->seen_generation=pool->generation;
        size_t *new_frozen=pool->new_frozen,*total_frozen=pool->total_frozen;
        float *min_raw_mu=pool->min_raw_mu;
        pthread_mutex_unlock(&pool->mutex);

        ExposureWorkspace *ws=&exp_workspace[arg->worker_id];
        memset(ws->partial_freq,0,nc*sizeof(fftwf_complex));
        ws->time_forward_ifft=ws->time_detector=ws->time_scatter_fft=0.0;
        for(int k=arg->worker_id;k<nexp;k+=n_exp_workspace)
            fused_exposure(k,ws,&new_frozen[k],&total_frozen[k],&min_raw_mu[k]);

        pthread_mutex_lock(&pool->mutex);
        if(--pool->pending==0)pthread_cond_signal(&pool->done_cv);
        pthread_mutex_unlock(&pool->mutex);
    }
    return NULL;
}
#endif

static void init_iteration_pool(void) {
#ifdef UPDC_PARALLEL
    if(n_exp_workspace>1) {
        memset(&iter_pool,0,sizeof(iter_pool));
        if(pthread_mutex_init(&iter_pool.mutex,NULL)||pthread_cond_init(&iter_pool.work_cv,NULL)||
           pthread_cond_init(&iter_pool.done_cv,NULL))fail("cannot initialize persistent iteration pool");
        iter_pool.threads=alloc((size_t)n_exp_workspace*sizeof(*iter_pool.threads));
        iter_pool.args=alloc((size_t)n_exp_workspace*sizeof(*iter_pool.args));
        for(int w=0;w<n_exp_workspace;++w) {
            iter_pool.args[w]=(IterWorkerArg){w,0,&iter_pool};
            if(pthread_create(&iter_pool.threads[w],NULL,persistent_iteration_worker,&iter_pool.args[w]))
                fail("cannot create persistent iteration worker");
        }
    }
#endif
}
static void destroy_iteration_pool(void) {
#ifdef UPDC_PARALLEL
    if(n_exp_workspace>1) {
        pthread_mutex_lock(&iter_pool.mutex);iter_pool.stop=1;pthread_cond_broadcast(&iter_pool.work_cv);pthread_mutex_unlock(&iter_pool.mutex);
        for(int w=0;w<n_exp_workspace;++w)if(pthread_join(iter_pool.threads[w],NULL))fail("cannot join persistent iteration worker");
        free(iter_pool.args);free(iter_pool.threads);
        pthread_cond_destroy(&iter_pool.work_cv);pthread_cond_destroy(&iter_pool.done_cv);pthread_mutex_destroy(&iter_pool.mutex);
        memset(&iter_pool,0,sizeof(iter_pool));
    }
#endif
}

static void run_fused_exposures(size_t *new_frozen,size_t *total_frozen,float *min_raw_mu) {
#ifdef UPDC_PARALLEL
    if(n_exp_workspace>1) {
        pthread_mutex_lock(&iter_pool.mutex);
        iter_pool.new_frozen=new_frozen;iter_pool.total_frozen=total_frozen;iter_pool.min_raw_mu=min_raw_mu;
        iter_pool.pending=n_exp_workspace;++iter_pool.generation;
        pthread_cond_broadcast(&iter_pool.work_cv);
        while(iter_pool.pending>0)pthread_cond_wait(&iter_pool.done_cv,&iter_pool.mutex);
        pthread_mutex_unlock(&iter_pool.mutex);return;
    }
#endif
    ExposureWorkspace *ws=&exp_workspace[0];memset(ws->partial_freq,0,nc*sizeof(fftwf_complex));
    ws->time_forward_ifft=ws->time_detector=ws->time_scatter_fft=0.0;
    for(int k=0;k<nexp;++k)fused_exposure(k,ws,&new_frozen[k],&total_frozen[k],&min_raw_mu[k]);
}

static void critical_worker_times(double *forward_ifft,double *detector,double *scatter_fft) {
    *forward_ifft=*detector=*scatter_fft=0.0;
    for(int w=0;w<n_exp_workspace;++w) {
        if(exp_workspace[w].time_forward_ifft>*forward_ifft)*forward_ifft=exp_workspace[w].time_forward_ifft;
        if(exp_workspace[w].time_detector>*detector)*detector=exp_workspace[w].time_detector;
        if(exp_workspace[w].time_scatter_fft>*scatter_fft)*scatter_fft=exp_workspace[w].time_scatter_fft;
    }
}

static void reduce_adjoint_frequency_to_correction(void) {
    for(size_t p=0;p<nc;++p) {
        float re=0.0f,im=0.0f;
        for(int w=0;w<n_exp_workspace;++w){re+=exp_workspace[w].partial_freq[p][0];im+=exp_workspace[w].partial_freq[p][1];}
        fft_freq[p][0]=re;fft_freq[p][1]=im;
    }
    inverse_crop(correction); /* one inverse FFT for sum_k H_k^T scatter_k */
}

void UPDC_ALL(const char *manifest,int width,int height,int input_width,int input_height,
              int exposures,int transpose_psf,int iterations,
              const char *iter_prefix,const char *weight_output,
              int fft_worker_count,int exposure_worker_count,int preprocess_worker_count,int plan_measure,int autotune,int retune_fft,
              const char *wisdom_path,int save_every_arg) {
    if(width<1||height<1||input_width<1||input_height<1||exposures<1||iterations<1||
       fft_worker_count<1||exposure_worker_count<1||preprocess_worker_count<1||save_every_arg<0)
        fail("invalid dimensions/exposures/iterations/workers/save_every");
    double total_start=wall_seconds(), stage_start, fft_autotune_time, fft_plan_time,
           psf_time, geometry_time, mapping_time, prewarm_time, data_cache_time, weight_time,
           init_time, sensitivity_time;
    double *iteration_times=alloc((size_t)iterations*sizeof(double));
    double *iter_estimate_fft=alloc((size_t)iterations*sizeof(double));
    double *iter_fused_wall=alloc((size_t)iterations*sizeof(double));
    double *iter_forward_ifft=alloc((size_t)iterations*sizeof(double));
    double *iter_fused_detector=alloc((size_t)iterations*sizeof(double));
    double *iter_scatter_fft=alloc((size_t)iterations*sizeof(double));
    double *iter_adjoint_reduce=alloc((size_t)iterations*sizeof(double));
    double *iter_update=alloc((size_t)iterations*sizeof(double));
    double *iter_output=alloc((size_t)iterations*sizeof(double));
    W=width;H=height;iw=input_width;ih=input_height;nexp=exposures;
    psf_transposed=transpose_psf?1:0;
    fft_workers=fft_worker_count;exposure_workers=exposure_worker_count;preprocess_workers=preprocess_worker_count;
    fft_plan_measure=plan_measure?1:0;fft_autotune_enabled=autotune?1:0;fft_retune=retune_fft?1:0;save_every=save_every_arg;
    fftw_wisdom_path=(wisdom_path&&wisdom_path[0])?wisdom_path:NULL;
    atomic_init(&fft_forward_count,0);atomic_init(&fft_inverse_count,0);
    if(!host_is_little_endian())fail("big-endian hosts are not supported by the raw float32 interchange format");
    if((size_t)W>SIZE_MAX/(size_t)H||(size_t)iw>SIZE_MAX/(size_t)ih)fail("image size overflow");
    nout=(size_t)W*H;nin=(size_t)iw*ih;
    if(nout>UINT32_MAX)fail("target grid exceeds 32-bit link index capacity");
    persistent_cache_enabled=(getenv("UPDC_PERSIST_CACHE") && strcmp(getenv("UPDC_PERSIST_CACHE"),"0")!=0);
    set_cache_root_from_manifest(manifest);
    printf("Persistent C preprocessing cache: %s\n",persistent_cache_enabled?"enabled":"disabled");fflush(stdout);
    parse_manifest(manifest);

    configure_fft_threads();
    if(fftw_wisdom_path) {
        int loaded=fftwf_import_wisdom_from_filename(fftw_wisdom_path);
        printf("FFTW wisdom %s: %s\n",fftw_wisdom_path,loaded?"loaded":"not found/unused");fflush(stdout);
    }

    int max_psf_width=0,max_psf_height=0;
    for(int k=0;k<nexp;++k) {
        if(expv[k].psf_width>max_psf_width)max_psf_width=expv[k].psf_width;
        if(expv[k].psf_height>max_psf_height)max_psf_height=expv[k].psf_height;
    }
    const int pad_w=psf_transposed?max_psf_height:max_psf_width;
    const int pad_h=psf_transposed?max_psf_width:max_psf_height;
    if(W>INT_MAX-pad_w+1||H>INT_MAX-pad_h+1)fail("FFT padding size overflow");
    int min_fw=W+pad_w-1;
    int min_fh=H+pad_h-1;
    stage_start=wall_seconds();
    select_fft_shape(min_fw,min_fh);
    fft_autotune_time=wall_seconds()-stage_start;
    nf=(size_t)FW*FH; nc=(size_t)FH*(FW/2+1);
    printf("Background-aware UPDC: %d exposures, input %dx%d, target %dx%d, FFT %dx%d, iterations %d, fft-workers %d, exposure-workers %d, preprocess-workers %d\n",
           nexp,iw,ih,W,H,FW,FH,iterations,fft_workers,exposure_workers,preprocess_workers); fflush(stdout);
    printf("Precision: f64 polygon coordinates/intersections; f32 area + numerical solver\n"); fflush(stdout);

    estimate=zalloc(nout*sizeof(float)); correction=zalloc(nout*sizeof(float));
    if((size_t)nexp>SIZE_MAX/nin)fail("freeze mask size overflow");
    mu_frozen=zalloc((size_t)nexp*nin*sizeof(uint8_t));
    denominator=zalloc(nout*sizeof(float)); scatter=zalloc(nout*sizeof(float));
    coverage=zalloc(nout*sizeof(float)); blur=alloc(nout*sizeof(float));
    fft_image=fftwf_alloc_real(nf); fft_freq=fftwf_alloc_complex(nc);
    estimate_freq=fftwf_alloc_complex(nc);
    if(!fft_image||!fft_freq||!estimate_freq)fail("FFT allocation failed");

    stage_start=wall_seconds();
    unsigned plan_flags=fft_plan_measure?FFTW_MEASURE:FFTW_ESTIMATE;
    forward_plan=fftwf_plan_dft_r2c_2d(FH,FW,fft_image,fft_freq,plan_flags);
    inverse_plan=fftwf_plan_dft_c2r_2d(FH,FW,fft_freq,fft_image,plan_flags);
    if(!forward_plan||!inverse_plan)fail("FFT planning failed");
    fft_plan_time=wall_seconds()-stage_start;
    printf("FFTW plan mode %s: %.3f s\n",fft_plan_measure?"MEASURE":"ESTIMATE",fft_plan_time);fflush(stdout);
    if(fftw_wisdom_path && !fftwf_export_wisdom_to_filename(fftw_wisdom_path))
        fprintf(stderr,"WARNING: could not write FFTW wisdom: %s\n",fftw_wisdom_path);

    stage_start=wall_seconds();
    prepare_all_psfs(max_psf_width,max_psf_height);
    psf_time=wall_seconds()-stage_start;
    printf("PSF spectra prepared in %.3f s (%d preprocess worker%s x %d FFT thread%s)\n",
           psf_time,preprocess_workers,preprocess_workers==1?"":"s",fft_workers,fft_workers==1?"":"s");fflush(stdout);

    stage_start=wall_seconds();
    build_area_links();
    geometry_time=wall_seconds()-stage_start;
    printf("Polygon link construction: %.3f s (%d exposure worker%s)\n",
           geometry_time,preprocess_workers,preprocess_workers==1?"":"s");fflush(stdout);

    stage_start=wall_seconds();
    size_t mapped_bytes=0;
    for(int k=0;k<nexp;++k) {
        linkmap[k]=open_links(k);
        mapped_bytes+=linkmap[k].bytes;
    }
    mapping_time=wall_seconds()-stage_start;
    printf("%d polygon LinkMaps persist through all iterations (%.3f GiB mapped)\n",
           nexp,(double)mapped_bytes/1073741824.0);fflush(stdout);

    stage_start=wall_seconds();load_input_planes();data_cache_time=wall_seconds()-stage_start;
    stage_start=wall_seconds();determine_weight_scale();weight_time=wall_seconds()-stage_start;
    stage_start=wall_seconds();SNA_HRA_stream();init_time=wall_seconds()-stage_start;
    stage_start=wall_seconds();precompute_denominator();sensitivity_time=wall_seconds()-stage_start;

    /* Allocate and first-touch all exposure workspaces before the mmap prewarm so
     * the final pre-iteration memory state is representative of the hot loop. */
    init_exposure_workspaces();
    prewarm_time=prewarm_mmaps();
    init_iteration_pool();

    size_t *new_by_exp=zalloc((size_t)nexp*sizeof(size_t));
    size_t *total_by_exp=zalloc((size_t)nexp*sizeof(size_t));
    float *min_by_exp=alloc((size_t)nexp*sizeof(float));

    char freeze_path[4096];
    if(snprintf(freeze_path,sizeof(freeze_path),"%s_mufreeze.csv",iter_prefix)>=(int)sizeof(freeze_path))
        fail("freeze CSV path is too long");
    FILE *freeze_csv=fopen(freeze_path,"w");if(!freeze_csv)fail("cannot write freeze CSV");
    fprintf(freeze_csv,"iteration,exposure,n_new_frozen,n_total_frozen,min_raw_mu,mu_marker,frozen_ratio\n");

    for(int iteration=1;iteration<=iterations;++iteration) {
        const double iter_start=wall_seconds();
        double t0=wall_seconds();
        prepare_estimate_spectrum();
        iter_estimate_fft[iteration-1]=wall_seconds()-t0;

        memset(new_by_exp,0,(size_t)nexp*sizeof(size_t));
        memset(total_by_exp,0,(size_t)nexp*sizeof(size_t));
        for(int k=0;k<nexp;++k)min_by_exp[k]=INFINITY;
        t0=wall_seconds();
        run_fused_exposures(new_by_exp,total_by_exp,min_by_exp);
        iter_fused_wall[iteration-1]=wall_seconds()-t0;
        critical_worker_times(&iter_forward_ifft[iteration-1],&iter_fused_detector[iteration-1],
                              &iter_scatter_fft[iteration-1]);

        size_t new_frozen=0,total_frozen=0;float min_raw_mu=INFINITY;
        for(int k=0;k<nexp;++k) {
            new_frozen+=new_by_exp[k];total_frozen+=total_by_exp[k];
            if(min_by_exp[k]<min_raw_mu)min_raw_mu=min_by_exp[k];
            fprintf(freeze_csv,"%d,%d,%zu,%zu,%.9g,%.9g,1\n",
                    iteration,k+1,new_by_exp[k],total_by_exp[k],(double)min_by_exp[k],(double)EPS_MU);
        }

        t0=wall_seconds();
        reduce_adjoint_frequency_to_correction();
        iter_adjoint_reduce[iteration-1]=wall_seconds()-t0;

        t0=wall_seconds();
        float sum=0.0f;size_t nnegative=0;
        for(size_t t=0;t<nout;++t) {
            if(denominator[t]>1e-8f) {
                float factor=correction[t]/denominator[t];
                if(!isfinite(factor))fail("non-finite multiplicative correction");
                if(factor<0) {
                    if(correction[t] < -1.0e-6f*max_sensitivity) {
                        fprintf(stderr,"negative correction at target %zu: C=%.9g numerator=%.9g denominator=%.9g\n",
                                t,(double)factor,(double)correction[t],(double)denominator[t]);
                        fail("substantially negative multiplicative correction");
                    }
                    factor=0.0f;
                }
                estimate[t]*=factor;
            } else estimate[t]=0;
            if(!isfinite(estimate[t]))fail("non-finite signed UPDC estimate");
            sum+=estimate[t];
            if(estimate[t]<0)nnegative++;
        }
        iter_update[iteration-1]=wall_seconds()-t0;

        t0=wall_seconds();
        const int should_save=(iteration==iterations)||(save_every>0 && iteration%save_every==0);
        if(should_save) {
            char iter_path[4096];
            if(snprintf(iter_path,sizeof(iter_path),"%s_%02d.bin",iter_prefix,iteration)>=(int)sizeof(iter_path))
                fail("iteration output path is too long");
            FILE *iter_file=fopen(iter_path,"wb");
            if(!iter_file)fail("cannot create iteration output");
            if(fwrite(estimate,sizeof(float),nout,iter_file)!=nout)fail("iteration output write failed");
            if(fclose(iter_file))fail("iteration output close failed");
            printf("UPDC iteration %d raw output: %s\n",iteration,iter_path);fflush(stdout);
        }
        iter_output[iteration-1]=wall_seconds()-t0;
        iteration_times[iteration-1]=wall_seconds()-iter_start;
        printf("UPDC iteration %d/%d: signed sum %.9g, negative %zu/%zu, frozen +%zu total %zu, min raw mu %.9g, wall %.3f s "
               "[estimateFFT %.3f fusedWall %.3f forwardIFFTcrit %.3f detectorCrit %.3f scatterFFTcrit %.3f adjointReduceIFFT %.3f update %.3f output %.3f]%s\n",
               iteration,iterations,(double)sum,nnegative,nout,new_frozen,total_frozen,(double)min_raw_mu,
               iteration_times[iteration-1],iter_estimate_fft[iteration-1],iter_fused_wall[iteration-1],iter_forward_ifft[iteration-1],
               iter_fused_detector[iteration-1],iter_scatter_fft[iteration-1],iter_adjoint_reduce[iteration-1],
               iter_update[iteration-1],iter_output[iteration-1],should_save?" SAVE":"");fflush(stdout);
    }
    if(fclose(freeze_csv))fail("freeze CSV close failed");

    for(size_t t=0;t<nout;++t)coverage[t]*=weight_scale;
    FILE *f=fopen(weight_output,"wb"); if(!f)fail("cannot create WHT output");
    if(fwrite(coverage,sizeof(float),nout,f)!=nout)fail("WHT output write failed");
    if(fclose(f))fail("WHT output close failed");

    char timing_path[4096];
    if(snprintf(timing_path,sizeof(timing_path),"%s_timing.csv",iter_prefix)>=(int)sizeof(timing_path))
        fail("timing output path is too long");
    FILE *timing=fopen(timing_path,"w");if(!timing)fail("cannot write timing CSV");
    fprintf(timing,"stage,wall_seconds\nfft_autotune,%.6f\nfft_plan,%.6f\npsf_cache,%.6f\n"
                   "polygon_links,%.6f\npersistent_mmap,%.6f\nmmap_prewarm,%.6f\ndata_cache,%.6f\nweight_scale,%.6f\n"
                   "initial_drizzle,%.6f\nsensitivity,%.6f\n",
            fft_autotune_time,fft_plan_time,psf_time,geometry_time,mapping_time,prewarm_time,
            data_cache_time,weight_time,init_time,sensitivity_time);
    for(int i=0;i<iterations;++i) {
        fprintf(timing,"iteration_%02d_estimate_fft,%.6f\n",i+1,iter_estimate_fft[i]);
        fprintf(timing,"iteration_%02d_fused_exposure_wall,%.6f\n",i+1,iter_fused_wall[i]);
        fprintf(timing,"iteration_%02d_forward_ifft_critical,%.6f\n",i+1,iter_forward_ifft[i]);
        fprintf(timing,"iteration_%02d_fused_detector_critical,%.6f\n",i+1,iter_fused_detector[i]);
        fprintf(timing,"iteration_%02d_scatter_fft_critical,%.6f\n",i+1,iter_scatter_fft[i]);
        fprintf(timing,"iteration_%02d_adjoint_reduce_ifft,%.6f\n",i+1,iter_adjoint_reduce[i]);
        fprintf(timing,"iteration_%02d_update,%.6f\n",i+1,iter_update[i]);
        fprintf(timing,"iteration_%02d_output,%.6f\n",i+1,iter_output[i]);
        fprintf(timing,"iteration_%02d_total,%.6f\n",i+1,iteration_times[i]);
    }
    fprintf(timing,"total_solver,%.6f\n",wall_seconds()-total_start);
    if(fclose(timing))fail("timing CSV close failed");

    const uint64_t nforward=atomic_load_explicit(&fft_forward_count,memory_order_relaxed);
    const uint64_t ninverse=atomic_load_explicit(&fft_inverse_count,memory_order_relaxed);
    printf("FFT transforms executed (excluding autotune/planning): forward=%llu inverse=%llu total=%llu\n",
           (unsigned long long)nforward,(unsigned long long)ninverse,
           (unsigned long long)(nforward+ninverse));fflush(stdout);

    struct rusage usage;
    if(!getrusage(RUSAGE_SELF,&usage)) {
#ifdef __APPLE__
        double peak_gib=(double)usage.ru_maxrss/1073741824.0;
#else
        double peak_gib=(double)usage.ru_maxrss/1048576.0;
#endif
        printf("Peak process RSS: %.3f GiB\n",peak_gib);fflush(stdout);
    }

    free(new_by_exp);free(total_by_exp);free(min_by_exp);
    destroy_iteration_pool();
    free_exposure_workspaces();
    for(int k=0;k<nexp;++k)close_links(&linkmap[k]);
    for(int k=0;k<nexp;++k) {
        if(expv[k].spectrum_map)munmap((void *)expv[k].spectrum_map,expv[k].spectrum_bytes);
        if(expv[k].spectrum)fclose(expv[k].spectrum);
    }
    fftwf_destroy_plan(forward_plan);fftwf_destroy_plan(inverse_plan);
    fftwf_free(fft_image);fftwf_free(fft_freq);fftwf_free(estimate_freq);
#ifdef UPDC_FFTW_THREADS
    if(fft_workers>1)fftwf_cleanup_threads();
#endif
    free(estimate);free(correction);free(denominator);free(scatter);free(coverage);
    free(blur);free(science_all);free(weight_all);free(mu_frozen);
    free(iteration_times);free(iter_estimate_fft);free(iter_fused_wall);free(iter_forward_ifft);free(iter_fused_detector);
    free(iter_scatter_fft);free(iter_adjoint_reduce);free(iter_update);free(iter_output);free(expv);free(linkmap);
}

static int self_test(void) {
    const UpdcPoint unit[4]={{1000.0,2000.0},{1000.0,2001.0},{1001.0,2001.0},{1001.0,2000.0}};
    const float area=updc_polygon_area(unit,4);
    if(fabsf(area-1.0f)>2.0e-6f) {
        fprintf(stderr,"self-test polygon area failed: %.9g\n",(double)area);
        return 1;
    }
    const UpdcPoint half[4]={{-0.5,-0.5},{-0.5,0.5},{0.0,0.5},{0.0,-0.5}};
    const float overlap=updc_overlap_area(half,0,0);
    if(fabsf(overlap-0.5f)>2.0e-6f) {
        fprintf(stderr,"self-test overlap failed: %.9g\n",(double)overlap);
        return 1;
    }
    if(sizeof(LinkMapHeader)!=64 || sizeof(TempLink)!=8) {
        fprintf(stderr,"self-test LinkMap ABI failed: header=%zu temp=%zu\n",sizeof(LinkMapHeader),sizeof(TempLink));
        return 1;
    }
    printf("UPDC C core self-test: OK\n");
    return 0;
}

int main(int argc,char **argv) {
    if(argc==2 && strcmp(argv[1],"--self-test")==0)return self_test();
    if(argc!=19)fail("usage: updc_jwst_c manifest width height input_width input_height exposures psf_transposed iterations iter_prefix wht.bin fft_workers exposure_workers preprocess_workers plan_measure autotune retune_fft wisdom_path save_every");
    UPDC_ALL(argv[1],atoi(argv[2]),atoi(argv[3]),atoi(argv[4]),atoi(argv[5]),
             atoi(argv[6]),atoi(argv[7]),atoi(argv[8]),argv[9],argv[10],
             atoi(argv[11]),atoi(argv[12]),atoi(argv[13]),atoi(argv[14]),atoi(argv[15]),atoi(argv[16]),argv[17],atoi(argv[18]));
    return 0;
}
