
/// \file FloatingPointComparison.hpp
/// \author jlippuner
/// \since Jan 9, 2014
///
/// \brief Defines functions to compare floating point numbers with
/// some tolerance.
///
///

#ifndef SKYNET_UTILITIES_FLOATINGPOINTCOMPARISON_HPP_
#define SKYNET_UTILITIES_FLOATINGPOINTCOMPARISON_HPP_

#include <math.h>

// copied from gsl/sys/fcmp.c so that the gsl library is not required to compile
int skynet_gsl_fcmp(const double x1, const double x2, const double epsilon);

/// function to compare doubles, returns:
///  0 if x = y (within tolerance)
/// +1 if x > y
/// -1 if x < y
int fcmp(const double x, const double y, const double epsilon = 1e-14);

#endif // SKYNET_UTILITIES_FLOATINGPOINTCOMPARISON_HPP_
