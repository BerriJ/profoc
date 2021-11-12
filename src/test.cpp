// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppClock)]]

#include <RcppArmadillo.h>
#include <RcppClock.h>
#include <thread>
using namespace arma;

// [[Rcpp::export]]
int test1()
{

    Rcpp::Clock clock;

    for (unsigned int x = 0; x < 100; x++)
    {
        mat A(1, 99);
        A.randn();
        mat B(99, 99);
        B.randn();
        mat C;
        clock.tick("1x99");
        C = A * B;
        clock.tock("1x99");
    }

    for (unsigned int x = 0; x < 100; x++)
    {
        mat A(2, 99);
        A.randn();
        mat B(99, 99);
        B.randn();
        mat C;
        clock.tick("2x99");
        C = A * B;
        clock.tock("2x99");
    }
    clock.stop("times");
    return 0;
}