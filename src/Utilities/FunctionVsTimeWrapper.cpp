/// \file FunctionVsTimeWrapper.cpp
/// \author jlippuner
/// \since Aug 31, 2014
///
/// \brief
///
///

#include "Utilities/FunctionVsTimeWrapper.hpp"

FunctionVsTimeWrapper::FunctionVsTimeWrapper(
    const std::function<double(double)>& function):
    mFunction(function) {}

double FunctionVsTimeWrapper::operator ()(const double time) const {
  return mFunction(time);
}

std::unique_ptr<FunctionVsTime<double>>
FunctionVsTimeWrapper::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(
      new FunctionVsTimeWrapper(*this));
}

std::unique_ptr<FunctionVsTime<double>> MakeFvT(
    const std::function<double(double)>& function) {
  return FunctionVsTimeWrapper(function).MakeUniquePtr();
}
