/// \file Homologous.hpp
/// \author jlippuner
/// \since Oct 6, 2017
///
/// \brief
///
///

#ifndef SKYNET_DENSITYPROFILES_HOMOLOGOUS_HPP_
#define SKYNET_DENSITYPROFILES_HOMOLOGOUS_HPP_

#include <memory>

#include "Utilities/FunctionVsTime.hpp"

class Homologous : public FunctionVsTime<double> {
public:
  Homologous(const double rho0, const double t0);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  double mRho0; // in g / cm^3
  double mT0; // in seconds
};

#endif // SKYNET_DENSITYPROFILES_HOMOLOGOUS_HPP_
