#ifndef profoc_h
#define profoc_h

#define RCPP_ARMADILLO_FIX_Field 1

#define OMP_PROC_BIND true
// #define OMP_DYNAMIC false

#ifdef _OPENMP
#include <omp.h>
#endif

#endif