/// \file PiecewiseLinearColorScale.cpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#include "MovieMaker/PiecewiseLinearColorScale.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

PieceWiseLinearColorScale::PieceWiseLinearColorScale(const double xMin,
    const double xMax, const bool useLog, const std::vector<double>& tics,
    const std::vector<std::string>& ticLabels,
    const std::vector<std::pair<double, double>>& red,
    const std::vector<std::pair<double, double>>& green,
    const std::vector<std::pair<double, double>>& blue,
    const std::vector<std::pair<double, double>>& alpha) :
    ColorScale(xMin, xMax, useLog, tics, ticLabels),
    mRed(red, false),
    mGreen(green, false),
    mBlue(blue, false),
    mAlpha(alpha, false) {

  assert(mRed.MinTime() == 0.0);
  assert(mRed.MaxTime() == 1.0);
  assert(mGreen.MinTime() == 0.0);
  assert(mGreen.MaxTime() == 1.0);
  assert(mBlue.MinTime() == 0.0);
  assert(mBlue.MaxTime() == 1.0);
  assert(mAlpha.MinTime() == 0.0);
  assert(mAlpha.MaxTime() == 1.0);
}

std::array<double, 4> PieceWiseLinearColorScale::GetRGBARange01(
    const double x) const {
  return { {mRed(x), mGreen(x), mBlue(x), mAlpha(x)}};
}

