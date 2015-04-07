/// \file FunctionVsTimeWrapper.hpp
/// \author jlippuner
/// \since Aug 31, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_FUNCTIONVSTIMEWRAPPER_HPP_
#define SKYNET_UTILITIES_FUNCTIONVSTIMEWRAPPER_HPP_

#include <functional>
#include <memory>

#include "Utilities/FunctionVsTime.hpp"

class FunctionVsTimeWrapper : public FunctionVsTime<double> {
public:
  FunctionVsTimeWrapper(const std::function<double(double)>& function);

  double operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<double>> MakeUniquePtr() const;

private:
  std::function<double(double)> mFunction;
};

std::unique_ptr<FunctionVsTime<double>> MakeFvT(
    const std::function<double(double)>& function);

#define FVT(func) MakeFvT(func).get()

#endif // SKYNET_UTILITIES_FUNCTIONVSTIMEWRAPPER_HPP_
