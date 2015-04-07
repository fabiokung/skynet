/// \file BisectionMethod.hpp
/// \author jonas
/// \since Apr 8, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_BISECTIONMETHOD_HPP_
#define SKYNET_UTILITIES_BISECTIONMETHOD_HPP_

#include <functional>

struct BisectionReturnValue {
  BisectionReturnValue(const double value, const bool isBestPossibleValue) :
      Value(value),
      IsBestPossibleValue(isBestPossibleValue) {}

  double Value;
  bool IsBestPossibleValue;
};

class BisectionMethod {
public:
  /// Return x such that f(x) = 0. x is assumed to be non-negative and f(x) is
  /// assumed to be monotonically increasing. xGuess is used to find the initial
  /// bracketing x values. xMin and xMax are the bounds of the valid range for
  /// x.
  static BisectionReturnValue FindRoot(
      const std::function<double(const double)>& f, const double xGuess,
      const double xMin, const double xMax, const int maxIter = 100,
      const double tolerance = 1.0E-14);

private:
  static BisectionReturnValue DoBisection(
      const std::function<double(const double)>& f, const double xLow,
      const double xHigh, const double fLow, const double fHigh,
      const int nSubintervals, const int maxIter, const double tolerance,
      const int iter);
};

#endif // SKYNET_UTILITIES_BISECTIONMETHOD_HPP_
