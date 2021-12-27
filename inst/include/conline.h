#ifndef conline_h
#define conline_h

#include <RcppArmadillo.h>

class conline
{
public:
    arma::mat y;
    arma::field<arma::cube> experts{0};
    const unsigned int &T = y.n_rows;
    const unsigned int &D = experts(0).n_rows;
    const unsigned int &P = experts(0).n_cols;
    const unsigned int &K = experts(0).n_slices;

    double x;

    // Constructors
    conline() = default; // Default constructor

    // Getters
    inline int getT() { return T; }
    inline int getD() { return D; }
    inline int getP() { return P; }
    inline int getK() { return K; }
};

#endif