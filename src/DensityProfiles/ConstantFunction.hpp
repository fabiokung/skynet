/// \file ConstantFunction.hpp
/// \author jlippuner
/// \since Oct 31, 2014
///
/// \brief
///
///

#ifndef SKYNET_DENSITYPROFILES_CONSTANTFUNCTION_HPP_
#define SKYNET_DENSITYPROFILES_CONSTANTFUNCTION_HPP_

#include <memory>

#include "Utilities/FunctionVsTime.hpp"

class ConstantFunction : public FunctionVsTime<double> {
public:
  ConstantFunction(const double value);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  double mValue;
};

#endif // SKYNET_DENSITYPROFILES_CONSTANTFUNCTION_HPP_
