/// \file Exp10.hpp
/// \author jlippuner
/// \since Nov 7, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_EXP10_HPP_
#define SKYNET_UTILITIES_EXP10_HPP_

#include <cmath>

inline double exp10(const double x) {
  return exp(log(10.0) * x);
}


#endif // SKYNET_UTILITIES_EXP10_HPP_
