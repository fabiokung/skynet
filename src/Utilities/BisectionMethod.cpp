/// \file BisectionMethod.cpp
/// \author jonas
/// \since Apr 8, 2014
///
/// \brief
///
///

#include "Utilities/BisectionMethod.hpp"

#include <stdexcept>
#include <vector>

#include "Utilities/CodeError.hpp"
#include "Utilities/FloatingPointComparison.hpp"

namespace { // unnamed so this can only be used in this file

union DoubleInt64 {
  double AsDouble;
  int64_t AsInt64;
};

}

//TODO make a Bisection exception

BisectionReturnValue BisectionMethod::FindRoot(
    const std::function<double(const double)>& f, const double xGuess,
    const double xMin, const double xMax, const int maxIter,
    const double tolerance) {
  double multiplier = 1.5;

  ASSERT(xMin < xMax, "xMin must be strictly less than xMax.");
  ASSERT(fcmp(f(xMin), 0.0, tolerance) <= 0, "f(xMin) must be <= 0");
  ASSERT(fcmp(f(xMax), 0.0, tolerance) >= 0, "f(xMax) must be >= 0");

  // check whether we already have a root
  double fGuess = f(xGuess);
  if (fcmp(fabs(fGuess), 0.0, tolerance) == 0)
    return BisectionReturnValue(xGuess, false);

  // find bracketing x values
  double xLow = 0.0;
  double xHigh = 0.0;

  if (fcmp(fGuess, 0.0, tolerance) > 0) {
    // guess is too large
    xHigh = xGuess;
    xLow = xGuess / multiplier;
    if (xLow < xMin)
        xLow = xMin;

    while (f(xLow) > 0.0) {
      if (xLow == xMin)
        throw std::runtime_error("f(xMin) = " + std::to_string(f(xMin))
            +" > 0, aborting Bisection");

      xLow /= multiplier;

      if (xLow < xMin)
        xLow = xMin;
    }
  } else {
    // guess is too small
    xHigh = xGuess * multiplier;
    if (xHigh > xMax)
      xHigh = xMax;
    xLow = xGuess;

    while (f(xHigh) < 0.0) {
      if (xHigh == xMax)
        throw std::runtime_error("f(xMax) = " + std::to_string(f(xHigh))
            + " < 0, aborting Bisection");

      xHigh *= multiplier;

      if (xHigh > xMax)
        xHigh = xMax;
    }
  }

  if (xHigh < xLow) {
    printf("xLow = %.20E, xHigh = %.20E\n", xLow, xHigh);
    throw std::runtime_error("xHigh < xLow");
  }

  return DoBisection(f, xLow, xHigh, f(xLow), f(xHigh), 10, maxIter, tolerance,
      0);
}

BisectionReturnValue BisectionMethod::DoBisection(
    const std::function<double(const double)>& f,
    const double xLow, const double xHigh, const double fLow,
    const double fHigh, const int nSubintervals, const int maxIter,
    const double tolerance, const int iter) {
  if (xHigh < xLow) {
    printf("xLow = %.20E, xHigh = %.20E\n", xLow, xHigh);
    throw std::runtime_error("xHigh < xLow");
  }

  //printf("iter = %i, fLow = %.20E, fHigh = %.20E\n", iter, fLow, fHigh);

  if (iter >= maxIter) {
    printf("# WARNING: Bisection did not converge\n");
    if (fabs(fLow) < fabs(fHigh))
      return BisectionReturnValue(xLow, false);
    else
      return BisectionReturnValue(xHigh, false);
  }

  DoubleInt64 xHighDoubleInt64, xLowDoubleInt64;
  xHighDoubleInt64.AsDouble = xHigh;
  xLowDoubleInt64.AsDouble = xLow;

  int64_t ulpDiff = xHighDoubleInt64.AsInt64 - xLowDoubleInt64.AsInt64;
  if (ulpDiff <= nSubintervals) {
    // we cannot split this range up into n subintervals, so evaluate f at all
    // possible values and pick the best
//    printf("ulp diff = %li\n, xLow = %.20E, xHigh = %.20E\n", ulpDiff, xLow,
//        xHigh);

    std::vector<double> xs(ulpDiff + 1);
    std::vector<double> fs(ulpDiff + 1);

    xs[0] = xLow;
    fs[0] = fLow;
    xs[ulpDiff] = xHigh;
    fs[ulpDiff] = fHigh;

    for (int i = 1; i < ulpDiff; ++i) {
      DoubleInt64 x;
      x.AsInt64 = xLowDoubleInt64.AsInt64 + i;
      xs[i] = x.AsDouble;
      fs[i] = f(xs[i]);
    }

    double xBest = xHigh;
    double fBest = fHigh;

    for (int i = 0; i < ulpDiff; ++i) {
      if (fabs(fs[i]) < fabs(fBest)) {
        fBest = fs[i];
        xBest = xs[i];
      }
    }

    //printf("num iter = %i\n", iter);
    return BisectionReturnValue(xBest, true);
  }

  // create sub intervals
  double dx = (xHigh - xLow) / (double)nSubintervals;
  std::vector<double> xs(nSubintervals + 1);
  std::vector<double> fs(nSubintervals + 1);

  xs[0] = xLow;
  fs[0] = fLow;
  xs[nSubintervals] = xHigh;
  fs[nSubintervals] = fHigh;

  for (int i = 1; i < nSubintervals; ++i) {
    xs[i] = xLow + dx * (double)i;
    if (xs[i] <= xs[i-1])
      throw std::runtime_error("xs are not monotonically increasing");
    fs[i] = f(xs[i]);
  }

  if (xs[nSubintervals] <= xs[nSubintervals - 1])
    throw std::runtime_error("xs are not monotonically increasing");

//  for (int i = 0; i <= nSubintervals; ++i) {
//    printf("x[%02i] = %30.20E, f = %30.20E\n", i, xs[i], fs[i]);
//  }

  // find new sub interval

  // check whether we found a solution
  for (int i = 0; i <= nSubintervals; ++i) {
    if (fabs(fs[i]) < tolerance) {
      //printf("num iter = %i\n", iter);
      return BisectionReturnValue(xs[i], false);
    }
  }

  double bestFAtSubintervalBoundary = fabs(fs[0]);
  int newIntervalLowerIdx = 0;

  for (int i = 1; i <= nSubintervals; ++i) {
    if (fs[i-1] * fs[i] <= 0.0) {
      // this is a subinterval containing a root
      double bestBoundaryF = std::min(fabs(fs[i-1]), fabs(fs[i]));
      if (bestBoundaryF < bestFAtSubintervalBoundary) {
        bestFAtSubintervalBoundary = bestBoundaryF;
        newIntervalLowerIdx = i - 1;
      }
    }
  }

  return DoBisection(f, xs[newIntervalLowerIdx], xs[newIntervalLowerIdx + 1],
      fs[newIntervalLowerIdx], fs[newIntervalLowerIdx + 1], nSubintervals,
      maxIter, tolerance, iter + 1);
}

