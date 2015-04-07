/// \file PowerLawContinuation.hpp
/// \author jlippuner
/// \since Sep 26, 2014
///
/// \brief
///
///

#ifndef SKYNET_DENSITYPROFILES_POWERLAWCONTINUATION_HPP_
#define SKYNET_DENSITYPROFILES_POWERLAWCONTINUATION_HPP_

#include <memory>

#include "Utilities/FunctionVsTime.hpp"
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class PowerLawContinuation : public FunctionVsTime<double> {
public:
  PowerLawContinuation(const PiecewiseLinearFunction& initialData,
      const double powerLawSlope, const double blendTimespan);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  // this piecewise linear function contains the initial data plus the data
  // region where the initial data and the power law are blended together
  PiecewiseLinearFunction mInitialData;

  // the raw power law starts at this time
  double mPowerLawStartTime;

  // these are the coefficients of the power law
  double mPowerLawSlope;
  double mPowerLawOffset;
};

#endif // SKYNET_DENSITYPROFILES_POWERLAWCONTINUATION_HPP_
