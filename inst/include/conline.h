#ifndef conline_h
#define conline_h

#include <RcppArmadillo.h>
#include "clock.h"

class conline
{
public:
    // Data
    arma::mat y;
    arma::field<arma::cube> experts;
    arma::vec tau;

    // Hyperparameters
    unsigned int lead_time = 0;
    int start;

    std::string loss_function = "quantile";
    double loss_parameter = 1.0;
    bool loss_gradient = true;

    std::string method = "bewa";

    std::unordered_map<std::string, arma::colvec> params;
    std::unordered_map<std::string, arma::colvec> params_basis_pr;
    std::unordered_map<std::string, arma::colvec> params_basis_mv;
    std::unordered_map<std::string, arma::colvec> params_hat_pr;
    std::unordered_map<std::string, arma::colvec> params_hat_mv;
    double forget_past_performance = 0.0;
    bool allow_quantile_crossing = false;
    bool trace = true;

    // Smoothing matrices
    arma::field<arma::sp_dmat> basis_pr;
    arma::field<arma::sp_dmat> basis_mv;
    arma::field<arma::sp_dmat> hat_pr;
    arma::field<arma::sp_dmat> hat_mv;

    // Starting vals
    arma::cube w0;
    arma::cube R0;

    // Externally specified loss/regret
    arma::field<arma::cube> loss_array;
    arma::field<arma::cube> regret_array;

    // Dimension shortcuts - for convenience
#define T y.n_rows
#define D experts(0).n_rows
#define P experts(0).n_cols
#define K experts(0).n_slices
#define X params.begin()->second.n_elem
#define T_E_Y int(experts.n_rows - T)

    // Internal objects

    // Loss
    arma::cube loss_for;
    arma::field<arma::cube> loss_exp;

    // Weights
    arma::field<arma::cube> weights_tmp;
    arma::field<arma::cube> weights;

    // Predictions
    arma::field<arma::cube> predictions_tmp;
    arma::cube predictions;

    // Performance related
    arma::mat chosen_params;
    arma::vec opt_index;
    arma::field<arma::cube> past_performance;
    arma::vec tmp_performance;
    arma::vec cum_performance;

    // Learning parameters
    arma::field<arma::cube> V;
    arma::field<arma::cube> E;
    arma::field<arma::cube> k;
    arma::field<arma::cube> eta;
    arma::field<arma::cube> R;
    arma::field<arma::cube> R_reg;
    arma::field<arma::cube> beta;
    arma::field<arma::cube> beta0field;

    // For benchmarking
    Rcpp::Clock clock;

    conline() = default; // Default constructor

    // Getters
    inline int getT() { return T; }
    inline int getD() { return D; }
    inline int getP() { return P; }
    inline int getK() { return K; }
    inline int getX() { return X; }
    // Methods
    void set_defaults();
    void set_grid_objects();
    void learn();
    void init_update(
        Rcpp::List &object,
        arma::mat &new_y,
        arma::field<arma::cube> &new_experts,
        const bool trace);
    Rcpp::List output();
    void teardown()
    {
        clock.stop();
    };
    ~conline() = default; // Default destructor
};

#endif