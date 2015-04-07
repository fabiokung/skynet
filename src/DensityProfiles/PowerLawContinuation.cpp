/// \file PowerLawContinuation.cpp
/// \author jlippuner
/// \since Sep 26, 2014
///
/// \brief
///
///

#include "DensityProfiles/PowerLawContinuation.hpp"

#include <cmath>
#include <stdexcept>

namespace { // unnamed so this function can only be used in this file

PiecewiseLinearFunction MakeInitialData(
    const PiecewiseLinearFunction& initialData, const double powerLawSlope,
    const double blendTimespan, double * pPowerLawStartTime,
    double * pPowerLawOffset) {
  // check the blend timespan is within the data range of the initial data
  double blendRegionStart = initialData.MaxTime() - blendTimespan;
  if (blendRegionStart < initialData.MinTime())
    throw std::invalid_argument("Tried to create a power law continuation "
        "with a blend time span that is bigger than initial data time span");

  *pPowerLawStartTime = initialData.MaxTime();

  // find power law offset by fitting straight line with slope powerLawSlope to
  // the data (in log-log space) in the blend region
  auto times = initialData.Times();
  auto values = initialData.Values();
  double sum = 0.0;
  int numPointsInBlendRegion = 0;

  for (unsigned int i = 0; i < times.size(); ++i) {
    if (times[i] < blendRegionStart)
      continue;

    sum += log(values[i]) - powerLawSlope * log(times[i]);
    ++numPointsInBlendRegion;
  }

  *pPowerLawOffset = sum / (double)numPointsInBlendRegion;

  // now do the blending
  for (unsigned int i = 0; i < times.size(); ++i) {
    if (times[i] < blendRegionStart)
      continue;

    // 0 at the beginning of the blend region, 1 at the end
    double frac = (times[i] - blendRegionStart) / blendTimespan;
    values[i] = (1.0 - frac) * values[i]
        + frac * exp(powerLawSlope * log(times[i]) + *pPowerLawOffset);
  }

  return PiecewiseLinearFunction(times, values, true);
}

}

PowerLawContinuation::PowerLawContinuation(
    const PiecewiseLinearFunction& initialData, const double powerLawSlope,
    const double blendTimespan) :
    mInitialData(MakeInitialData(initialData, powerLawSlope, blendTimespan,
        &mPowerLawStartTime, &mPowerLawOffset)),
    mPowerLawSlope(powerLawSlope) {}

double PowerLawContinuation::operator()(const double time) const {
  if (time < mPowerLawStartTime)
    return mInitialData(time);
  else
    return exp(mPowerLawSlope * log(time) + mPowerLawOffset);
}

std::unique_ptr<FunctionVsTime<double>>
PowerLawContinuation::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(
      new PowerLawContinuation(*this));
}
