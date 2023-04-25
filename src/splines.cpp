#include <splines.h>

#include <misc.h>

#include <RcppArmadillo.h>
// include header file from splines2 package
#include <splines2Armadillo.h>

//' @title Create B-Spline basis
//'
//' @description This function creates a B-Spline matrix.
//'
//' @param x Vector of values.
//' @param knots Vector of knots.
//' @param deg Degree of the Spline functions.
//' @param periodic Whether the basis should be periodic or not.
//' @param intercept Whether the firs column should be kept.
//' @return Returns a matrix of B-Spline basis functions.
//' @examples
//' n <- 9
//' deg <- 3
//' mu <- 0.35
//' x <- 0:1000 / 1000
//'
//' knots <- make_knots(n, mu = mu, deg = deg)
//'
//' B <- splines2_basis(x, knots, deg)
//' ts.plot(B, col = 1:dim(B)[2])
//'
//' # Periodic Case
//' B <- splines2_basis(x, knots, deg, periodic = TRUE)
//' ts.plot(B, col = 1:dim(B)[2])
//'
//' @export
// [[Rcpp::export]]
arma::mat splines2_basis(const arma::vec &x,
                         const arma::vec &knots,
                         const unsigned int deg,
                         const bool &periodic = false,
                         const bool &intercept = true)
{
    arma::mat B; // Basis matrix

    if (!periodic)
    {
        splines2::BSpline bs_obj{x, deg, knots};
        B = bs_obj.basis(true);
    }
    else
    {
        // We will only use the inner and boundary knots for the periodic case
        unsigned int order = deg + 1;
        arma::uvec inner_idx = arma::regspace<arma::uvec>(order,
                                                          knots.n_elem - order - 1);
        arma::uvec bound_idx = {deg, knots.n_elem - order};

        splines2::PeriodicBSpline bs_obj{x, knots(inner_idx), deg, knots(bound_idx)};
        B = bs_obj.basis(true);
    }

    if (!intercept)
        B.shed_col(0);

    return B;
}

using namespace arma;

//' @title B-Spline penalty
//'
//' @description This function calculates the B-Spline basis penalty.
//' It follows the procedure outlined in the paper by Zheyuan Li, Jiguo
//' Cao, 2022 "General P-Splines for Non-Uniform B-Splines"
//' \doi{10.48550/arXiv.2201.06808}.
//' For equidistant knots it coincides with the usual penalty based
//' on the identitiy. For non-equidistant knots it is a weighted penalty
//' with respect to the knot distances.
//' In addition to the above, we added the possibility to calculate
//' periodic penalties which are based on the periodic differencing matrices.
//'
//' @param knots Vector of knots.
//' @param order Order of the Basis (degree + 1).
//' @param periodic Whether the penalties should be periodic or not.
//' @param max_diff Maximum difference order to calculate.
//'
//' @return Returns a list of (order - 1) penalty matrices.
//'
//' @examples
//' \dontrun{
//' # Equidistant knots with order 2
//' knots <- 1:10
//'
//' P <- penalty(knots, order = 2)
//'
//' print(P[[1]]) # First differences
//'
//' # Non equidistant knots
//' knots <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
//'
//' P <- penalty(knots, order = 4)
//'
//' print(P[[1]]) # First differences
//' print(P[[2]]) # Second differences
//' print(P[[3]]) # Third differences
//'
//' # Periodic penalty for equidistant knots
//' oder <- 4
//' deg <- order - 1
//' knots <- 1:15
//'
//' penalty(knots, order = order, periodic = TRUE)[[1]]
//' penalty(knots, order = order, periodic = TRUE)[[2]]
//' penalty(knots, order = order, periodic = TRUE)[[3]]
//' }
//'
//' @export
// [[Rcpp::export]]
arma::field<arma::sp_mat> penalty(
    const arma::vec &knots,
    const unsigned int &order,
    const bool &periodic = false,
    const unsigned int &max_diff = 999)
{
    int K = knots.n_elem;
    unsigned int i = 1;
    arma::field<arma::sp_mat> D(order);
    arma::field<arma::sp_mat> P(std::min(order - 1, max_diff));

    if (!periodic)
    {

        D(0) = eye(K - order, K - order);
        mat d = diff(eye(K - order, K - order));

        // While i < order calculate the next difference matrix and the respective scaled penalty
        while (i < order && i <= max_diff)
        {
            arma::vec h = diff_cpp(knots.rows(i, K - 1 - i), order - i, 1) / (order - i);
            mat W_inv = diagmat(1 / h);
            D(i) = W_inv * d.submat(0, 0, d.n_rows - i, d.n_cols - i) * D(i - 1);
            P(i - 1) = D(i).t() * D(i);
            P(i - 1) *= std::pow(arma::mean(h), 2 * i);
            i++;
        }
    }
    else
    {
        if (K - 2 * order < 3)
        {
            throw std::invalid_argument("At least order-1 inner knots are needed for periodic penalties. \n K <- length(knots) \n J <- K - 2*order \n J >= 3 # Must be true");
        }

        D(0) = eye(K - 2 * order + 1, K - 2 * order + 1);
        mat dp = diff(eye(K - 2 * order + 2, K - 2 * order + 2));
        dp.shed_col(0);
        dp(0, dp.n_cols - 1) = -1;

        // Extend the knot sequence
        arma::uvec inner_idx = arma::regspace<arma::uvec>(order,
                                                          K - order - 1);
        arma::uvec bound_idx = {order - 1, K - order};

        // We need this sequence to calculate the weights
        arma::vec knots_ext = knots.subvec(bound_idx(0), bound_idx(1));

        knots_ext = join_cols(knots_ext,
                              knots(inner_idx.head(order - 1)) + knots(bound_idx(1)) - knots(bound_idx(0)));

        K = knots_ext.n_elem; // Update number of Knots

        arma::vec h = diff_cpp(knots_ext.rows(0, K - 2), order - 1, 1);
        h /= (order - 1);
        arma::mat w_inv = diagmat(1 / h);

        while (i < order && i <= max_diff)
        {
            D(i) = w_inv * dp * D(i - 1);
            P(i - 1) = D(i).t() * D(i);
            P(i - 1) *= std::pow(arma::mean(h), 2 * i);
            i++;
        }
    }

    return P;
}

// [[Rcpp::export]]
arma::mat periodic_adjacency(const int &size)
{
    arma::mat adj(size, size, arma::fill::zeros);

    for (int i = 0; i < size; ++i)
    {
        adj(i, (i + 1) % size) = 1;
        adj(i, (i - 1 + size) % size) = 1;
    }

    return adj;
}

// [[Rcpp::export]]
arma::mat adjacency_to_incidence(const arma::mat &adj)
{
    int cols = adj.n_cols;

    int edge = 0;
    arma::mat incidence(0, cols);

    for (int col = 0; col < cols; ++col)
    {
        // We only look at half the adjacency matrix, so that we only add each
        // edge to the incidence matrix once.
        for (int row = col; row >= 0; --row)
        {
            if (adj(row, col) > 0)
            {
                incidence.resize(cols, edge + 1);
                incidence(row, edge) = 1;
                incidence(col, edge) = 1;
                ++edge;
            }
        }
    }

    return incidence;
}

// [[Rcpp::export]]
arma::sp_mat make_basis_matrix(const arma::vec &x,
                               const arma::vec &knots,
                               const unsigned int deg,
                               const bool &periodic = false)
{
    mat B;

    if (knots.n_elem == 1)
    {
        mat B_(x.n_elem, 1, fill::ones);
        B = B_;
    }
    else
    {
        B = splines2_basis(x, knots, deg, periodic, true);
    }

    B.clean(1E-10);

    // Remove columns without contribution
    B = B.cols(find(sum(B) >= 1E-6));

    return sp_mat(B);
}

// [[Rcpp::export]]
arma::sp_mat make_hat_matrix(
    const arma::vec &x,
    const arma::vec &knots,
    const int deg,
    const double &bdiff,
    const double &lambda,
    const bool &periodic)
{

    mat B; // Basis matrix
    mat P; // Penalty matrix

    B = splines2_basis(x, knots, deg, periodic, true);

    // Field of penalties up to second differences
    arma::field<arma::sp_mat> Ps = penalty(knots, deg + 1, periodic, 2);

    // It may look strange that we check for deg > 1 here, but penalties for
    // higher differences can only be computed if deg > 1.
    if (deg > 1)
    {
        P = (2 - bdiff) * Ps(0) + (bdiff - 1) * Ps(1);
    }
    else
    {
        P = Ps(0);
    }

    mat H = B * arma::pinv(B.t() * B + lambda * P) * B.t(); // Hat matrix
    H.clean(1E-10);

    return sp_mat(H);
}
