/// \file StopWatch.cpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#include "Utilities/StopWatch.hpp"

#include <stdexcept>

StopWatch::StopWatch() :
    mElapsedTime(duration_t::zero()),
    mIsRunning(false),
    mStartTime() {
}

void StopWatch::Start() {
  auto time = clock_t::now();

  if (mIsRunning)
    throw std::logic_error("Tried to start running stop watch.");

  mIsRunning = true;
  mStartTime = time;
}

void StopWatch::Stop() {
  auto time = clock_t::now();

  if (!mIsRunning)
    throw std::logic_error("Tried to stop clock that is not running.");

  mIsRunning = false;
  duration_t elapsedTime = std::chrono::duration_cast<duration_t>(
      time - mStartTime);
  mElapsedTime += elapsedTime;
}

void StopWatch::Reset() {
  if (mIsRunning)
    throw std::logic_error("Tried to reset running stop watch.");
  mElapsedTime = duration_t::zero();
}

double StopWatch::GetElapsedTimeInSeconds() const {
  typedef std::chrono::duration<double, std::ratio<1, 1>> sec_t;
  return std::chrono::duration_cast<sec_t>(mElapsedTime).count();
}

double StopWatch::GetElapsedTimeInMilliseconds() const {
  typedef std::chrono::duration<double, std::ratio<1, 1000>> millisec_t;
  return std::chrono::duration_cast<millisec_t>(mElapsedTime).count();
}

double StopWatch::GetElapsedTimeInMicroseconds() const {
  typedef std::chrono::duration<double, std::ratio<1, 1000000>> microsec_t;
  return std::chrono::duration_cast<microsec_t>(mElapsedTime).count();
}

double StopWatch::GetElapsedTimeInNanoseconds() const {
  typedef std::chrono::duration<double, std::ratio<1, 1000000000>> nanosec_t;
  return std::chrono::duration_cast<nanosec_t>(mElapsedTime).count();
}
