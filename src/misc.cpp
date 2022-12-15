#include <RcppArmadillo.h>
#include <random>
#include <set>

// [[Rcpp::export]]
std::set<unsigned long long int> sample_int(
    unsigned long long int N, 
    unsigned long long int size)
{
    std::mt19937 rng(std::random_device{}());

    // Create an empty set of integers.
    std::set<unsigned long long int> set;

    while (set.size() < size)
    {
        unsigned long long int value = std::uniform_int_distribution<int>(1, N)(rng);
        set.insert(value);
    }

    return set;
}