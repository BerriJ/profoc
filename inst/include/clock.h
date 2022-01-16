#ifndef clock_h
#define clock_h

#include <RcppArmadillo.h>
#include <chrono>
#include <string>
#include <map>

#ifndef _OPENMP
inline int omp_get_thread_num() { return 0; }
#endif

namespace Rcpp
{
    namespace sc = std::chrono;

    class Clock
    {
        using tp = sc::high_resolution_clock::time_point;
        using keypair = std::pair<std::string, int>;
        using timesmap = std::map<keypair, tp>;

    private:
        std::string name;
        timesmap tickmap;
        std::vector<double> timers;
        std::vector<std::string> names;

    public:
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

        // Ppass data to R
        void stop()
        {
            DataFrame df = DataFrame::create(
                Named("Name") = names,
                Named("Nanoseconds") = timers);
            df.attr("class") = "cppclock";
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