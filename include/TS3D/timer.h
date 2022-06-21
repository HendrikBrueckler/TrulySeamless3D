#ifndef TIME_H
#define TIME_H

#include <chrono>
#include <vector>

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_start);
double subtractTimes(std::chrono::high_resolution_clock::time_point& t_end,
                     std::chrono::high_resolution_clock::time_point& t_start);

double getAvgTime(std::vector<double> time);
double getTotalTime(std::vector<double> time);

#endif // TIME_H
