/// \file ColorScale.cpp
/// \author jlippuner
/// \since Jun 12, 2014
///
/// \brief
///
///

#include "MovieMaker/ColorScale.hpp"

#include <cassert>

ColorScale::ColorScale(const double xMin, const double xMax, const bool useLog,
    const std::vector<double>& tics, const std::vector<std::string>& ticLabels):
    mScale(xMin, xMax, useLog, tics, ticLabels) {  }

std::array<double, 4> ColorScale::GetRGBA(const double x) const {
  return GetRGBARange01(mScale(x));
}
