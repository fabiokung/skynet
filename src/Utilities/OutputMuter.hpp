/// \file OutputMuter.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SRC_UTILITIES_OUTPUTMUTER_HPP_
#define SRC_UTILITIES_OUTPUTMUTER_HPP_

#include <cstdio>

class OutputMuter {
public:
  static void Mute();
  static void Unmute();

private:
  OutputMuter() { }

  static bool mIsCurrentlyMuted;
  static fpos_t mOldStdoutPos;
  static int mOldStdoutFileDescriptor;
};

#endif // SRC_UTILITIES_OUTPUTMUTER_HPP_
