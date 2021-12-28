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

        typedef std::map<std::pair<std::string, int>,
                         std::chrono::high_resolution_clock::time_point>
            TimesMap;

    private:
        std::vector<double> timers;
        size_t ticks, tocks = 0;
        TimesMap tickmap, tockmap;
        std::pair<std::string, int> key;

    public:
        // start a timer
        void tick(std::string name)
        {
            key.first = name;
            key.second = ticks;
            tickmap.insert(
                std::pair<std::pair<std::string, int>, std::chrono::high_resolution_clock::time_point>(key, std::chrono::high_resolution_clock::now()));
            ticks += 1;
        }

        // stop a timer
        void tock(std::string name)
        {
            key.first = name;
            key.second = tocks;
            tockmap.insert(
                std::pair<std::pair<std::string, int>, std::chrono::high_resolution_clock::time_point>(key, std::chrono::high_resolution_clock::now()));
            tocks += 1;
        }

        // calculate timer durations
        void stop(std::string var_name)
        {
            std::vector<std::string> keys;
            keys.reserve(tickmap.size());
            std::vector<std::chrono::high_resolution_clock::time_point> ticks;
            ticks.reserve(tickmap.size());
            for (auto kv : tickmap)
            {
                keys.push_back(kv.first.first);
                ticks.push_back(kv.second);
            }

            std::vector<std::chrono::high_resolution_clock::time_point> tocks;
            tocks.reserve(tockmap.size());
            for (auto kv : tockmap)
            {
                tocks.push_back(kv.second);
            }

            for (std::size_t i = 0; i < ticks.size(); ++i)
            {
                timers.push_back(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        tocks[i] - ticks[i])
                        .count());
            }
            DataFrame df = DataFrame::create(
                Named("Name") = keys,
                Named("Nanoseconds") = timers);
            Environment env = Environment::global_env();
            env[var_name] = df;
        }
    };
}

#endif