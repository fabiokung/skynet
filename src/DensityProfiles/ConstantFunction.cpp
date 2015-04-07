/// \file ConstantFunction.cpp
/// \author jlippuner
/// \since Oct 31, 2014
///
/// \brief
///
///

#include "DensityProfiles/ConstantFunction.hpp"

ConstantFunction::ConstantFunction(const double value):
    mValue(value) {}

double ConstantFunction::operator()(const double /*time*/) const {
  return mValue;
}

std::unique_ptr<FunctionVsTime<double>> ConstantFunction::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<double>>(new ConstantFunction(mValue));
}
