/// \file Homologous.cpp
/// \author jlippuner
/// \since Oct 6, 2017
///
/// \brief
///
///

#include "DensityProfiles/Homologous.hpp"

#include <cmath>
#include <stdexcept>

Homologous::Homologous(const double rho0, const double t0):
    mRho0(rho0),
    mT0(t0) {
  if (mRho0 <= 0.0)
    throw std::invalid_argument("Rho0 for an Homologous density profile must "
        "be positive");
  if (mT0 <= 0.0)
    throw std::invalid_argument("t0 for an Homologous density profile must "
        "be positive");
}

double Homologous::operator()(const double time) const {
  return mRho0 * pow(mT0 / time, 3.0);
}

std::unique_ptr<FunctionVsTime<double>> Homologous::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(new Homologous(mRho0, mT0));
}
