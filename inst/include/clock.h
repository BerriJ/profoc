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
        using keypair = std::pair<std::string, unsigned int>;
        using timesmap = std::map<keypair, tp>;

    private:
        timesmap tickmap;
        std::vector<double> averages;
        std::vector<unsigned int> counts;
        std::vector<std::string> unique_names;

    public:
        std::string name;
        std::vector<unsigned long long int> timers;
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

            // Loop over unique names
            for (unsigned int i = 0; i < unique_names.size(); i++)
            {
                unsigned long long int sum = 0;
                unsigned long int count = 0;

                // Loop over all names
                for (unsigned long int j = 0; j < names.size(); j++)
                {
                    if (names[j] == unique_names[i])
                    {
                        // Sum up all timers with the same name
                        sum += timers[j];
                        count++;
                    }
                }

                // Calculate average, round to 3 decimal places,
                // and convert from microseconds to milliseconds
                averages.push_back(std::round(sum / double(count)) / 1e+3);

                counts.push_back(count);
            }
        }

        // Pass data to R / Python
        void stop()
        {
            aggregate();

            // DataFrame df = DataFrame::create(
            //     Named("Name") = unique_names,
            //     Named("Milliseconds") = averages,
            //     Named("Count") = counts);
            // Environment env = Environment::global_env();
            // env[name] = df;
        }

        // Destructor
        ~Clock()
        {
            stop();
        }
    };
}
#endif