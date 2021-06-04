#ifndef pmax_arma_h
#define pmax_arma_h

#include <common_header.h>

mat pmax_arma(const mat &x, const double &bound);

mat pmin_arma(const mat &x, const double &bound);

vec diff_cpp(vec x, unsigned int lag, unsigned int differences);

mat get_combinations(const mat &x, const vec &y);

#endif