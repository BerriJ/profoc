#ifndef clock_h
#define clock_h

#include <Rcpp.h>
#include <chrono>
#include <string>
#include <vector>
#include <map>

#ifndef _OPENMP
inline int omp_get_thread_num() { return 0; }
#endif

namespace Rcpp
{
    template <typename T>
    void remove_duplicates(std::vector<T> &vec)
    {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }

    namespace sc = std::chrono;

    class Clock
    {
        using tp = sc::high_resolution_clock::time_point;
        using keypair = std::pair<std::string, unsigned int>;
        using timesmap = std::map<keypair, tp>;

    private:
        std::string name;                      // Name of R object to return
        timesmap tickmap;                      // Map of start times
        std::vector<std::string> names,        // Vector of identifiers
            unique_names;                      // Vector of unique identifiers
        std::vector<unsigned long int> counts; // Count occurence of identifiers
        std::vector<double> means, sds;        // Output vecs of mean and sd
        std::vector<unsigned long long int>    // Observed durations
            timers;

    public:
        // Init - Set name of R object
        Clock() : name("times") {}
        Clock(std::string name) : name(name) {}

        // start a timer - save time
        void tick(std::string &&name)
        {
            keypair key(std::move(name), omp_get_thread_num());

#pragma omp critical
            tickmap[key] = sc::high_resolution_clock::now();
        }

        // stop a timer - calculate time difference and save key
        void
        tock(std::string &&name)
        {
            keypair key(std::move(name), omp_get_thread_num());

#pragma omp critical
            {
                timers.push_back(
                    sc::duration_cast<sc::microseconds>(
                        sc::high_resolution_clock::now() -
                        tickmap[key])
                        .count());
                names.push_back(std::move(key.first));
            }
        }

        // Pass data to R / Python
        void aggregate()
        {
            // Create copy of names called unique_names
            unique_names = names;
            remove_duplicates(unique_names);

            for (unsigned int i = 0; i < unique_names.size(); i++)
            {
                unsigned long int count = 0;
                double mean = 0, M2 = 0, variance = 0;

                for (unsigned long int j = 0; j < names.size(); j++)
                {
                    if (names[j] == unique_names[i])
                    {
                        // Welford's online algorithm for mean and variance
                        double delta = timers[j] - mean;
                        count++;
                        mean += delta / count;
                        M2 += delta * (timers[j] - mean) * 1e-3;
                    }
                }

                // Save count
                counts.push_back(count);

                // Save average, round to 3 decimal places
                means.push_back(std::round(mean) * 1e-3);

                // Calculate sample variance
                variance = M2 / (count);
                // Save standard deviation, round to 3 decimal places
                sds.push_back(
                    std::round(std::sqrt(variance * 1e-3) * 1e+3) * 1e-3);
            }
        }

        // Pass data to R / Python
        void stop()
        {
            aggregate();

            DataFrame df = DataFrame::create(
                Named("Name") = unique_names,
                Named("Milliseconds") = means,
                Named("SD") = sds,
                Named("Count") = counts);
            Environment env = Environment::global_env();
            env[name] = df;
        }

        // Destructor
        ~Clock()
        {
            stop();
        }
    };
}
#endif