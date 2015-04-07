/// \file ExpTMinus3.cpp
/// \author jlippuner
/// \since Sep 4, 2014
///
/// \brief
///
///

#include "DensityProfiles/ExpTMinus3.hpp"

#include <cmath>
#include <stdexcept>

ExpTMinus3::ExpTMinus3(const double rho0, const double tau):
    mRho0(rho0),
    mTau(tau) {
  if (mRho0 <= 0.0)
    throw std::invalid_argument("Rho0 for an ExpTMinus3 denisty profile must "
        "be positive");
  if (mTau <= 0.0)
    throw std::invalid_argument("Tau for an ExpTMinus3 denisty profile must "
        "be positive");
}

double ExpTMinus3::operator()(const double time) const {
  if (time <= 3.0 * mTau)
    return mRho0 * exp(-time / mTau);
  else
    return mRho0 * pow(3.0 * mTau / (time * exp(1.0)), 3.0);
}

std::unique_ptr<FunctionVsTime<double>> ExpTMinus3::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(new ExpTMinus3(mRho0, mTau));
}
