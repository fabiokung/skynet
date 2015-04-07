/// \file FloatingPointExceptions.cpp
/// \author jlippuner
/// \since Aug 19, 2014
///
/// \brief
///
///

#include "Utilities/FloatingPointExceptions.hpp"

void FloatingPointExceptions::Enable() {
#ifdef SKYNET_ENABLE_FPE
  // first clear any existing exceptions
  feclearexcept(mExceptions);
  feenableexcept(mExceptions);
#endif // SKYNET_ENABLE_FPE
}

void FloatingPointExceptions::Disable() {
#ifdef SKYNET_ENABLE_FPE
  fedisableexcept(mExceptions);
#endif // SKYNET_ENABLE_FPE
}
