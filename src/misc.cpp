#include <Rcpp.h>
#include <cstdint>
#include <random>
#include <set>

// [[Rcpp::export]]
std::set<uint64_t> sample_int(
    uint64_t N,
    uint64_t size,
    uint32_t seed)
{
    std::mt19937_64 rng(seed);

    // Create an empty set of integers.
    std::set<uint64_t> set;

    while (set.size() < size)
    {
        uint64_t value = std::uniform_int_distribution<uint64_t>(1, N)(rng);
        set.insert(value);
    }

    return set;
}
