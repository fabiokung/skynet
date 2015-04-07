/// \file Scale.cpp
/// \author jonas
/// \since Jun 14, 2014
///
/// \brief
///
///

#include "MovieMaker/Scale.hpp"

#include <cassert>
#include <cmath>

Scale::Scale(const double xMin, const double xMax, const bool useLog,
    const std::vector<double>& tics,
    const std::vector<std::string>& ticLabels):
    mXMin(useLog ? log(xMin) : xMin),
    mXMax(useLog ? log(xMax) : xMax),
    mDx(mXMax - mXMin),
    mUseLog(useLog),
    mTics(tics),
    mTicLabels(ticLabels) {
  assert(xMin < xMax);
  assert((!useLog) || (xMin > 0.0));
  assert(tics.size() == ticLabels.size());
}

double Scale::operator()(const double x) const {
  double xIn = mXMin;
  if (mUseLog && (x > 0.0)) {
    xIn = log(x);
  }
  if (!mUseLog)
    xIn = x;

  if (xIn < mXMin)
    xIn = mXMin;
  if (xIn > mXMax)
    xIn = mXMax;

  return (xIn - mXMin) / mDx;
}
