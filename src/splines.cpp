#include <splines.h>

#include <misc.h>

#include <RcppArmadillo.h>
// include header file from splines2 package
#include <splines2Armadillo.h>

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
        if (!intercept)
            B.shed_col(0);
    }
    else
    {
        unsigned int order = deg + 1;
        arma::uvec inner_idx = arma::regspace<arma::uvec>(order,
                                                          knots.n_elem - order - 1);
        arma::uvec bound_idx = {deg, knots.n_elem - order};

        // Create periodic mspline object
        splines2::PeriodicMSpline ps(x, knots(inner_idx), deg, knots(bound_idx));
        B = ps.basis(true);

        if (!intercept)
            B.shed_col(0);

        // We will use https://doi.org/10.6339/21-JDS1020 eq. (1) to convert
        // the mspline basis to a bspline basis

        // We need this sequence to calculate the weights
        arma::vec knots_ext = knots.subvec(bound_idx(0), bound_idx(1));
        knots_ext = arma::join_cols(knots_ext,
                                    knots(inner_idx.head(deg)) + knots(bound_idx(1)));

        for (unsigned int i = 0; i < B.n_cols; i++)
        {
            double w = knots_ext(1 - intercept + i + order) -
                       knots_ext(1 - intercept + i);
            B.col(i) *= w / order;
        }
    }

    return B;
}

using namespace arma;

// [[Rcpp::export]]
arma::vec make_knots(const double &kstep, const double &a = 1, const int deg = 3, const bool &even = false)
{
    vec x;
    vec xa;
    vec xb;

    if (even)
    {
        x = arma::linspace((0 - deg * 2 * kstep), 1 - kstep, (1 / (2 * kstep) - 0.5) + deg + 1);
        xa = sign(x) % pow(arma::abs(x), a) / 2;
        xb = 1 - reverse(xa.subvec(0, xa.n_elem - 1));
    }
    else
    {
        x = arma::linspace(0 - deg * 2 * kstep, 1, 1 / (2 * kstep) + deg + 1);
        xa = sign(x) % pow(arma::abs(x), a) / 2;
        xb = 1 - reverse(xa.subvec(0, xa.n_elem - 2));
    }

    return join_cols(xa, xb);
}

//' @title B-Spline penalty
//'
//' @description This function calculates the B-Spline basis penalty.
//' It follows the procedure outlined in the paper by Zheyuan Li, Jiguo
//' Cao, 2022 "General P-Splines for Non-Uniform B-Splines"
//' \doi{10.48550/arXiv.2201.06808}.
//' For equidistant knots it coincides with the usual penalty based
//' on the identitiy. For non-equidistant knots it is a weighted penalty
//' with respect to the knot distances.
//'
//' @param knots Vector of knots.
//' @param order Order of the Basis (degree + 1).
//' @param max_diff Maximum difference order to calculate.
//'
//' @return Returns a list of (order - 1) penalty matrices.
//'
//' @examples
//' \dontrun{
//' # Equidisan knots with order 2
//' knots <- 1:10
//'
//' P <- penalty(knots, order = 2)
//'
//' print(P[[1]]) # First differences
//'
//' # Non-equidistant knots
//' knots <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
//'
//' P <- penalty(knots, order = 4)
//'
//' print(P[[1]]) # First differences
//' print(P[[2]]) # Second differences
//' print(P[[3]]) # Third differences
//' }
//'
//' @export
// [[Rcpp::export]]
arma::field<arma::sp_mat> penalty(
    const arma::vec &knots,
    const int &order,
    const int &max_diff = 999)
{
    int K = knots.n_elem;
    arma::field<arma::sp_mat> D(order);
    arma::field<arma::sp_mat> P(std::min(order - 1, max_diff));
    D(0) = eye(K - order, K - order);
    mat d = diff(eye(K - order, K - order));

    int i = 1;

    // While i < order calculate the next difference matrix and the respective scaled penalty
    while (i <= std::min(order - 1, max_diff))
    {
        arma::vec h = diff_cpp(knots.rows(i, K - 1 - i), order - i, 1) / (order - i);
        mat W_inv = diagmat(1 / h);
        D(i) = W_inv * d.submat(0, 0, d.n_rows - i, d.n_cols - i) * D(i - 1);
        P(i - 1) = D(i).t() * D(i);
        P(i - 1) *= std::pow(arma::mean(h), 2 * i);
        i++;
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
arma::mat penalty_periodic(
    const arma::vec &knots,
    const unsigned int &order)
{

    // TODO merge into penalty function

    unsigned int K = knots.n_elem;

    // Create incidence from adjacency matrix
    arma::mat adj = periodic_adjacency(K - 2 * order + 1);
    arma::mat inc = adjacency_to_incidence(adj);

    // Extend the knot sequence
    arma::uvec inner_idx = arma::regspace<arma::uvec>(order,
                                                      K - order - 1);
    arma::uvec bound_idx = {order - 1, K - order};

    // We need this sequence to calculate the weights
    arma::vec knots_ext = knots.subvec(bound_idx(0), bound_idx(1));
    knots_ext = join_cols(knots_ext,
                          knots(inner_idx.head(order - 1)) + knots(bound_idx(1)));

    K = knots_ext.n_elem; // Update number of Knots

    int i = 1;

    arma::vec h = diff_cpp(knots_ext.rows(i - 1, K - 1 - i), order - i, 1);
    h /= (order - i);
    arma::mat w_inc = diagmat(1 / h) * inc;
    arma::mat P = -w_inc.t() * w_inc;
    P.diag() *= -1;

    P *= std::pow(arma::mean(h), 2 * i);

    return P;
}

// [[Rcpp::export]]
arma::sp_mat make_hat_matrix(
    const arma::vec &x,
    const double &kstep,
    const double &lambda,
    const double &bdiff,
    const int deg,
    const double &a,
    const bool &even)
{
    mat H;

    if (kstep <= 0.5)
    {
        vec knots = make_knots(kstep, a, deg, even);
        int m = knots.n_elem - 2 * (deg)-2; // Number of inner knots

        mat B = splines2_basis(x, knots, deg, false, true);

        mat P(m + deg + 1, m + deg + 1);

        // Field of penalties up to second differences
        arma::field<arma::sp_mat> Ps = penalty(knots, deg + 1, 2);

        if (deg > 1)
        {
            P = (2 - bdiff) * Ps(0) + (bdiff - 1) * Ps(1);
        }
        else
        {
            P = Ps(0);
        }

        H = B * arma::pinv(B.t() * B + lambda * P) * B.t();
        H.clean(1E-10);
    }
    else
    {
        mat identity(x.n_elem, x.n_elem, fill::eye);
        H = identity;
    }

    arma::sp_mat sp_H = sp_mat(H); // Return hat matrix
    return sp_H;
}

// [[Rcpp::export]]
arma::sp_mat make_basis_matrix(const arma::vec &x, const double &kstep, const int deg, const double &a, const bool &even)
{
    mat B;

    // Will be passed to make_hat_matrix
    if (kstep <= 0.5)
    {
        vec knots = make_knots(kstep, a, deg, even);
        B = splines2_basis(x, knots, deg, false, true);
        // Remove columns without contribution
        B = B.cols(find(sum(B) >= 1E-6));
        B.clean(1E-10);
    }
    else
    {
        mat B_(x.n_elem, 1, fill::ones);
        B = B_;
    }

    sp_mat out(B);

    return out;
}

// [[Rcpp::export]]
arma::sp_mat make_basis_matrix2(const arma::vec &x,
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
arma::sp_mat make_hat_matrix2(
    const arma::vec &x,
    const arma::vec &knots,
    const int deg,
    const double &bdiff,
    const double &lambda,
    const bool &periodic)
{

    mat B; // Basis matrix
    mat P; // Penalty matrix

    if (!periodic)
    {
        B = splines2_basis(x, knots, deg, periodic, true);

        // Field of penalties up to second differences
        arma::field<arma::sp_mat> Ps = penalty(knots, deg + 1, 2);

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
    }
    else
    {
        if (bdiff != 1)
        {
            throw std::invalid_argument("Currently only bdiff 1 is supported for periodic splines.");
        }

        B = splines2_basis(x, knots, deg, periodic, true);
        P = penalty_periodic(knots, deg + 1);
    }
    mat H = B * arma::pinv(B.t() * B + lambda * P) * B.t(); // Hat matrix
    H.clean(1E-10);

    return sp_mat(H);
}
