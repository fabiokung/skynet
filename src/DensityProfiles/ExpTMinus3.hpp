/// \file ExpTMinus3.hpp
/// \author jlippuner
/// \since Sep 4, 2014
///
/// \brief
///
///

#ifndef SKYNET_DENSITYPROFILES_EXPTMINUS3_HPP_
#define SKYNET_DENSITYPROFILES_EXPTMINUS3_HPP_

#include <memory>

#include "Utilities/FunctionVsTime.hpp"

class ExpTMinus3 : public FunctionVsTime<double> {
public:
  ExpTMinus3(const double rho0, const double tau);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  double mRho0; // in g / cm^3
  double mTau; // in seconds
};

#endif // SKYNET_DENSITYPROFILES_EXPTMINUS3_HPP_
