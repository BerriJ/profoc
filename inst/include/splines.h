#ifndef splines_h
#define splines_h

#include <common_header.h>

vec make_knots(const double &kstep, const double &a = 1, const int deg = 3);

mat splineDesign_rcpp(const vec &x, const vec &knots, const int &deg);

mat make_difference_matrix(const vec &knots, const int &bdiff, const int deg);

mat make_hat_matrix(const vec &x,
                    const double &kstep,
                    const double &lambda,
                    const double &bdiff,
                    const int deg,
                    const double &a);

vec spline_fit(const vec &y,
               const vec &x,
               const double &lambda = 1,
               const int &ndiff = 1,
               const int &deg = 3,
               const double &knot_distance = 0.1,
               const double &knot_distance_power = 1);

#endif