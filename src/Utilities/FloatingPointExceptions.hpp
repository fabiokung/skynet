/// \file FloatingPointExceptions.hpp
/// \author jlippuner
/// \since Aug 19, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_FLOATINGPOINTEXCEPTIONS_HPP_
#define SKYNET_UTILITIES_FLOATINGPOINTEXCEPTIONS_HPP_

// find out whether we are using a GNU compiler, use information given at
// http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros

#if (defined(__GNUC__) || defined(__GNUG__)) && !defined(__APPLE__)
#define SKYNET_ENABLE_FPE
#else
#undef SKYNET_ENABLE_FPE
#warning "Floating point exceptions are only supported with GNU on Linux"
#endif

#ifdef SKYNET_ENABLE_FPE
#include <cfenv>
#endif // SKYNET_ENABLE_FPE

class FloatingPointExceptions {
public:
  static void Enable();

  static void Disable();

private:
  // we don't enable FE_UNDERFLOW, since that happens regularly when calcuting
  // rates that are very small
#ifdef SKYNET_ENABLE_FPE
  static const int mExceptions = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW;
#else
  static const int mExceptions = 0;
#endif // SKYNET_ENABLE_FPE
};

#endif // SKYNET_UTILITIES_FLOATINGPOINTEXCEPTIONS_HPP_
