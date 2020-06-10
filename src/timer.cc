#include "TrulySeamless3D/timer.h"

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_end,
                     std::chrono::high_resolution_clock::time_point& t_start)
{
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
    return elapsed.count() / 1000000.0;
}

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_start)
{
    auto t_end = std::chrono::high_resolution_clock::now();
    return subtractTimes(t_end, t_start);
}

double getAvgTime(std::vector<double> time)
{
    return getTotalTime(time) / time.size();
}

double getTotalTime(std::vector<double> time)
{
    double t = 0.0;
    for (unsigned int i = 0; i < time.size(); i++)
        t += time[i];
    return t;
}
