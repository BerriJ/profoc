#ifndef clock_h
#define clock_h

#include <RcppArmadillo.h>
#include <chrono>
#include <string>
#include <map>

namespace Rcpp
{
    class Clock
    {

    private:
        std::map<std::string,
                 std::chrono::high_resolution_clock::time_point>
            tickmap;
        std::string key;
        std::vector<std::string> keys;
        std::vector<double> timers;

    public:
        // start a timer - save time
        void tick(std::string name)
        {
            tickmap.insert(
                std::pair<std::string, std::chrono::high_resolution_clock::time_point>(name, std::chrono::high_resolution_clock::now()));
        }

        // stop a timer - calculate time difference and save key
        void tock(std::string name)
        {
            timers.push_back(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now() - tickmap[name])
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