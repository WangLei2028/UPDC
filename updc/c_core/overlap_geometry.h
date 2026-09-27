#ifndef UPDC_OVERLAP_GEOMETRY_H
#define UPDC_OVERLAP_GEOMETRY_H

/*
 * Geometry policy for the optimized C release:
 *   - polygon vertices and clipping/intersection coordinates stay in double;
 *   - polygon areas are accumulated in local float coordinates.
 *
 * Translating each polygon to its first vertex before casting to float avoids
 * catastrophic cancellation from ~5k-pixel absolute target coordinates while
 * keeping the area hot path single precision.
 */
typedef struct { double x, y; } UpdcPoint;

float updc_polygon_area(const UpdcPoint *polygon, int count);
float updc_overlap_area(const UpdcPoint quad[4], int x, int y);

#endif
