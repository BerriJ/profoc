#ifndef clock_h
#define clock_h

#include <RcppArmadillo.h>
#include <chrono>
#include <string>
#include <map>

namespace Rcpp
{
    namespace sc = std::chrono;

    class Clock
    {
        using tp = sc::high_resolution_clock::time_point;
        using timesmap = std::map<std::string, tp>;

    private:
        timesmap tickmap;
        std::string key;
        std::vector<std::string> keys;
        std::vector<double> timers;

    public:
        // start a timer - save time
        void tick(std::string name)
        {
            tickmap.insert(
                std::pair<std::string, tp>(name,
                                           sc::high_resolution_clock::now()));
        }

        // stop a timer - calculate time difference and save key
        void tock(std::string name)
        {
            timers.push_back(
                sc::duration_cast<sc::nanoseconds>(
                    sc::high_resolution_clock::now() - tickmap[name])
                    .count());
            keys.push_back(name);
        }

        // calculate timer durations
        void stop(std::string var_name)
        {
            DataFrame df = DataFrame::create(
                Named("Name") = keys,
                Named("Nanoseconds") = timers);
            df.attr("class") = "cppclock";
            Environment env = Environment::global_env();
            env[var_name] = df;
        }
    };
}
#endif