#ifndef clock_h
#define clock_h

#include <RcppArmadillo.h>
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
        using keypair = std::pair<std::string, int>;
        using timesmap = std::map<keypair, tp>;

    private:
        timesmap tickmap;

    public:
        std::string name;
        std::vector<double> timers;
        std::vector<std::string> names;

        // Init - Set name of R object
        Clock() : name("times") {}
        Clock(std::string name_) : name(name_) {}

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
                    sc::duration_cast<sc::nanoseconds>(
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
            std::vector<std::string> unique_names = names;
            remove_duplicates(unique_names);

            std::vector<std::tuple<std::string, double, int>>
                table(unique_names.size());

            std::vector<double> averages(unique_names.size());
            std::vector<int> counts(unique_names.size());

            // Loop over unique names
            for (unsigned int i = 0; i < unique_names.size(); i++)
            {
                int sum = 0;
                int count = 0;

                // Loop over all names
                for (unsigned int j = 0; j < names.size(); j++)
                {
                    if (names[j] == unique_names[i])
                    {
                        sum += timers[j];
                        count++;
                    }
                }

                // Calculate average, convert to milliseconds, round to 3 dec
                averages[i] = (std::round((sum * 1e-3) / double(count)) / 1e+3);
                counts[i] = count;
            }

            DataFrame df = DataFrame::create(
                Named("Name") = unique_names,
                Named("Milliseconds") = averages,
                Named("Count") = counts);
            Environment env = Environment::global_env();
            env[name] = df;
        }

        void stop()
        {
            aggregate();
        }

        // Destructor
        ~Clock()
        {
            stop();
        }
    };
}
#endif