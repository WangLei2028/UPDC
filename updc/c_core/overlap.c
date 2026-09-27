#include "overlap_geometry.h"

#include <math.h>

typedef UpdcPoint Point;

enum Boundary {
    BOUNDARY_LEFT = 0,
    BOUNDARY_RIGHT = 1,
    BOUNDARY_BOTTOM = 2,
    BOUNDARY_TOP = 3
};

static inline int inside(Point p, enum Boundary boundary, double value) {
    switch (boundary) {
        case BOUNDARY_LEFT:   return p.x >= value;
        case BOUNDARY_RIGHT:  return p.x <= value;
        case BOUNDARY_BOTTOM: return p.y >= value;
        case BOUNDARY_TOP:    return p.y <= value;
    }
    return 0;
}

static inline Point intersection(Point a, Point b, enum Boundary boundary, double value) {
    Point out = a;
    double delta;
    if (boundary == BOUNDARY_LEFT || boundary == BOUNDARY_RIGHT) {
        delta = b.x - a.x;
        if (fabs(delta) > 1.0e-15) {
            const double t = (value - a.x) / delta;
            out.x = value;
            out.y = a.y + t * (b.y - a.y);
        }
    } else {
        delta = b.y - a.y;
        if (fabs(delta) > 1.0e-15) {
            const double t = (value - a.y) / delta;
            out.x = a.x + t * (b.x - a.x);
            out.y = value;
        }
    }
    return out;
}

static int clip_boundary(
    const Point *input,
    int input_count,
    Point *output,
    enum Boundary boundary,
    double value
) {
    int output_count = 0;
    Point previous = input[input_count - 1];
    int previous_inside = inside(previous, boundary, value);
    for (int i = 0; i < input_count; ++i) {
        const Point current = input[i];
        const int current_inside = inside(current, boundary, value);
        if (current_inside) {
            if (!previous_inside) {
                output[output_count++] = intersection(previous, current, boundary, value);
            }
            output[output_count++] = current;
        } else if (previous_inside) {
            output[output_count++] = intersection(previous, current, boundary, value);
        }
        previous = current;
        previous_inside = current_inside;
    }
    return output_count;
}

float updc_polygon_area(const Point *polygon, int count) {
    /*
     * Clipping coordinates remain double.  Translate the polygon to a local
     * origin before casting to float so area products stay O(pixel^2), not
     * O(absolute_coordinate^2).  This keeps f32 area stable on ~5k mosaics.
     */
    if (count < 3) return 0.0f;
    const double origin_x = polygon[0].x;
    const double origin_y = polygon[0].y;
    float twice_area = 0.0f;
    for (int i = 0; i < count; ++i) {
        const int j = (i + 1) % count;
        const float xi = (float)(polygon[i].x - origin_x);
        const float yi = (float)(polygon[i].y - origin_y);
        const float xj = (float)(polygon[j].x - origin_x);
        const float yj = (float)(polygon[j].y - origin_y);
        twice_area += fmaf(xi, yj, -(xj * yi));
    }
    return 0.5f * fabsf(twice_area);
}

float updc_overlap_area(const Point quad[4], int x, int y) {
    Point first[16];
    Point second[16];
    Point *input = first;
    Point *output = second;
    int count = 4;
    const enum Boundary boundaries[4] = {
        BOUNDARY_LEFT, BOUNDARY_RIGHT, BOUNDARY_BOTTOM, BOUNDARY_TOP
    };
    const double values[4] = {
        (double)x - 0.5, (double)x + 0.5,
        (double)y - 0.5, (double)y + 0.5
    };

    for (int j = 0; j < 4; ++j) input[j] = quad[j];
    for (int i = 0; i < 4; ++i) {
        count = clip_boundary(input, count, output, boundaries[i], values[i]);
        if (count < 3) return 0.0f;
        Point *swap = input;
        input = output;
        output = swap;
    }
    return updc_polygon_area(input, count);
}
