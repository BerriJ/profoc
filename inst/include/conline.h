#ifndef conline_h
#define conline_h

#include <RcppArmadillo.h>

class conline
{
public:
    conline() {} // Default initializer, somehow important for RcppModules
    arma::mat A;
    double x;
    double y;
};

#endif