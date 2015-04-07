/// \file StopWatch.hpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_STOPWATCH_HPP_
#define SKYNET_UTILITIES_STOPWATCH_HPP_

#include <chrono>

class StopWatch {
public:
  typedef std::chrono::nanoseconds duration_t;
  typedef std::chrono::high_resolution_clock clock_t;

  StopWatch();

  void Start();

  void Stop();

  void Reset();

  double GetElapsedTimeInSeconds() const;

  double GetElapsedTimeInMilliseconds() const;

  double GetElapsedTimeInMicroseconds() const;

  double GetElapsedTimeInNanoseconds() const;

  bool IsRunning() const {
    return mIsRunning;
  }

private:
  duration_t mElapsedTime;
  bool mIsRunning;
  clock_t::time_point mStartTime;
};

#endif // SKYNET_UTILITIES_STOPWATCH_HPP_
