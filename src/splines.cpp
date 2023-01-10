#include <splines.h>

#include <misc.h>

#include <RcppArmadillo.h>
// include header file from splines2 package
#include <splines2Armadillo.h>

static arma::mat splines2_basis(const arma::vec &x,
                                const arma::vec &knots,
                                const unsigned int deg)
{
    splines2::BSpline bs_obj{x, deg, knots};
    return bs_obj.basis(true);
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

// Expose splines::splineDesign to Rcpp
// static mat splineDesign_rcpp(const vec &x, const vec &knots, const int &deg)
// {
//     Rcpp::Environment pkg = Rcpp::Environment::namespace_env("splines");
//     Rcpp::Function f = pkg["splineDesign"];
//     mat y = Rcpp::as<arma::mat>(f(knots, x, deg + 1, 0, true));
//     return y;
// }

// [[Rcpp::export]]
arma::mat wt_delta(const arma::vec &h)
{
    int r = h.n_elem;
    arma::vec pos_neg = {-1, 1};
    arma::vec x = repmat(pos_neg, r, 1) % repelem(1 / h, 2, 1);
    arma::uvec i = regspace<uvec>(0, r - 1);
    i = repelem(i, 2, 1);
    arma::uvec p(r + 2, fill::zeros);
    p.rows(1, p.n_rows - 2) = regspace<uvec>(1, 2, 2 * r - 1);
    p.row(p.n_rows - 1) = 2 * r;
    arma::mat D(r, r + 1);
    int col_ptr = 0;
    for (unsigned int iter = 0; iter < x.n_elem; iter++)
    {
        D(i(iter), col_ptr) = x(iter);
        if (p(col_ptr + 1) == iter + 1)
        {
            col_ptr++;
        }
    }
    return D;
}

//' @title Difference matrix for B-Spline penalty
//'
//' @description This function calculates the difference matrices for the
//' B-Spline basis penalty. It follows the procedure outlined in the paper
//' by Zheyuan Li, Jiguo Cao, 2022 "General P-Splines for Non-Uniform B-Splines"
//' \doi{10.48550/arXiv.2201.06808}.
//' For equidistant knots it coincides with the usual difference matrix based
//' on the identitiy. For non-equidistant knots it is a weighted difference
//' with respect to the knot distances.
//'
//' @param knots Vector of knots.
//' @param order Order of the Basis (degree + 1).
//' @param max_diff Maximum difference order to calculate.
//'
//' @return Returns a list of (order -1) difference matrices for
//' computing the B-Spline penalty.
//'
//' @examples
//' \dontrun{
//' # Equidisan knots with order 2
//' knots <- 1:10
//'
//' D <- penalty(knots, order = 2)
//'
//' print(D[[1]]) # First differences
//'
//' # Non-equidistant knots
//' knots <- c(0, 0, 0, 0, 1, 3, 4, 4, 4, 4)
//'
//' D1 <- penalty(knots, order = 4)
//'
//' print(D1[[1]]) # First differences
//' print(D1[[2]]) # Second differences
//' print(D1[[3]]) # Third differences
//' }
//'
//' @export
// [[Rcpp::export]]
arma::field<arma::mat> penalty(
    const arma::vec &knots,
    const int &order,
    const int &max_diff = 999)
{

    int i = 1;
    int K = knots.n_elem;
    arma::field<arma::mat> D(std::min(order - 1, max_diff));
    arma::field<arma::mat> P(D);
    arma::vec h = diff_cpp(knots.rows(i, K - 1 - i), order - i, 1) / (order - i);

    D(0) = wt_delta(h);
    P(0) = D(0).t() * D(0);
    P(0) *= std::pow(arma::mean(h), 2 * i);

    i++;

    // While i < order calculate the next difference matrix and save it into row i-1 of D and increment i
    while (i <= std::min(order - 1, max_diff))
    {
        h = diff_cpp(knots.rows(i, K - 1 - i), order - i, 1) / (order - i);
        D(i - 1) = wt_delta(h) * D(i - 2);
        P(i - 1) = D(i - 1).t() * D(i - 1);
        P(i - 1) *= std::pow(arma::mean(h), 2 * i);
        i++;
    }

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

        mat B = splines2_basis(x, knots, deg);

        mat P(m + deg + 1, m + deg + 1);

        // Field of penalties up to second differences
        arma::field<arma::mat> Ps = penalty(knots, deg + 1, 2);

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
        B = splines2_basis(x, knots, deg);
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
                                const unsigned int deg)
{
    mat B;

    if (knots.n_elem == 1)
    {
        mat B_(x.n_elem, 1, fill::ones);
        B = B_;
    }
    else
    {
        B = splines2_basis(x, knots, deg);
        // Remove columns without contribution
        B = B.cols(find(sum(B) >= 1E-6));
        B.clean(1E-10);
    }

    sp_mat out(B);

    return out;
}

// [[Rcpp::export]]
arma::sp_mat make_hat_matrix2(
    const arma::vec &x,
    const arma::vec &knots,
    const int deg,
    const double &bdiff,
    const double &lambda)
{

    mat H;

    int m = knots.n_elem - 2 * (deg)-2; // Number of inner knots

    mat B = splines2_basis(x, knots, deg);

    mat P(m + deg + 1, m + deg + 1);

    // Field of penalties up to second differences
    arma::field<arma::mat> Ps = penalty(knots, deg + 1, 2);

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

    arma::sp_mat sp_H = sp_mat(H); // Return hat matrix
    return sp_H;
}
