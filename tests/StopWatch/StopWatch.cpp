/// \file StopWatch.cpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#include "Utilities/StopWatch.hpp"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <thread>

#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  StopWatch sw;

  bool caughtException = false;
  try {
    sw.Stop();
  } catch (...) {
    caughtException = true;
  }
  if (!caughtException) {
    printf("Failed to catch stopping stopped stop watch.\n");
    return EXIT_FAILURE;
  }

  sw.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(350));
  sw.Stop();

  std::this_thread::sleep_for(std::chrono::milliseconds(500));

  sw.Start();
  std::this_thread::sleep_for(std::chrono::milliseconds(150));
  sw.Stop();

  printf("%.20e s\n", sw.GetElapsedTimeInSeconds());
  printf("%.20e ms\n", sw.GetElapsedTimeInMilliseconds());
  printf("%.20e micro s\n", sw.GetElapsedTimeInMicroseconds());
  printf("%.20e ns\n", sw.GetElapsedTimeInNanoseconds());

  if (fabs(0.5 - sw.GetElapsedTimeInSeconds()) / 0.5 > 0.01) {
    printf("Failed to measure time in seconds.\n");
    return EXIT_FAILURE;
  }

  if (fabs(500.0 - sw.GetElapsedTimeInMilliseconds()) / 500.0 > 0.01) {
    printf("Failed to measure time in milliseconds.\n");
    return EXIT_FAILURE;
  }

  if (fabs(500000.0 - sw.GetElapsedTimeInMicroseconds()) / 500000.0 > 0.01) {
    printf("Failed to measure time in microseconds.\n");
    return EXIT_FAILURE;
  }

  if (fabs(500000000.0 - sw.GetElapsedTimeInNanoseconds()) / 500000000.0
      > 0.01) {
    printf("Failed to measure time in nanoseconds.\n");
    return EXIT_FAILURE;
  }

  sw.Start();
  caughtException = false;
  try {
    sw.Start();
  } catch (...) {
    caughtException = true;
  }
  if (!caughtException) {
    printf("Failed to catch starting running stop watch.\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


