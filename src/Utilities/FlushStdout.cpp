/// \file FlushStdout.cpp
/// \author jlippuner
/// \since Aug 14, 2014
///
/// \brief
///
///

#include "Utilities/FlushStdout.hpp"

#include <chrono>
#include <cstdio>
#include <thread>

void FlushStdout::Flush() {
  // flush stdout
  fflush(stdout);

  // need to wait a little to make sure stdout really gets flushed
  std::this_thread::sleep_for(std::chrono::nanoseconds(1));
}
