/// \file FunctionVsTime.hpp
/// \author jlippuner
/// \since Aug 31, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_FUNCTIONVSTIME_HPP_
#define SKYNET_UTILITIES_FUNCTIONVSTIME_HPP_

#include <memory>

template<typename T>
class FunctionVsTime {
public:
  virtual ~FunctionVsTime() {}

  virtual T operator()(const double time) const =0;

  virtual std::unique_ptr<FunctionVsTime<T>> MakeUniquePtr() const =0;
};

#endif // SKYNET_UTILITIES_FUNCTIONVSTIME_HPP_
