#ifndef conline_h
#define conline_h

#include <RcppArmadillo.h>

class conline
{
public:
    // Data
    arma::mat y;
    arma::field<arma::cube> experts{2, 3, 4};
    arma::vec tau;

    // Hyperparameters
    unsigned int lead_time = 0;

    std::string loss_function{"quantile"};
    double loss_parameter;
    bool loss_gradient;

    std::string method;

    Rcpp::NumericMatrix param_grid;
    double forget_past_performance;
    bool allow_quantile_crossing;
    bool trace;

    // Smoothing matrices
    arma::field<arma::sp_dmat> basis_pr;
    arma::field<arma::sp_dmat> basis_mv;
    arma::field<arma::sp_dmat> hat_pr;
    arma::field<arma::sp_dmat> hat_mv;

    // Starting vals
    arma::cube w0;
    arma::cube R0;

    // Externally specified loss/regret
    arma::field<arma::cube> loss_array{0};
    arma::field<arma::cube> regret_array{0};

    // Dimension references - for convenience
    const unsigned int &T = y.n_rows;
    const unsigned int &D = experts(0).n_rows;
    const unsigned int &P = experts(0).n_cols;
    const unsigned int &K = experts(0).n_slices;
    const unsigned int &X = param_grid.rows();

    // Constructors
    conline() = default; // Default constructor

    // Getters
    inline int getT() { return T; }
    inline int getD() { return D; }
    inline int getP() { return P; }
    inline int getK() { return K; }
    inline int getX() { return X; }

    // Methods
    void init_objects();

    // Destructors
};

#endif