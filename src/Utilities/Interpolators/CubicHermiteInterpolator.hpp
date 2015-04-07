/// \file CubicHermiteInterpolator.hpp
/// \author jlippuner
/// \since Jan 7, 2015
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_CUBICHERMITEINTERPOLATOR_HPP_
#define SKYNET_UTILITIES_CUBICHERMITEINTERPOLATOR_HPP_

#include "Utilities/Interpolators/InterpolatorBase.hpp"

template<typename T>
class CubicHermiteInterpolator : public InterpolatorBase<T> {
public:
  CubicHermiteInterpolator(const std::vector<double>& times,
      const std::vector<T>& values);

  CubicHermiteInterpolator(
      const std::vector<std::pair<double, T>>& points);

  std::unique_ptr<FunctionVsTime<T>> MakeUniquePtr() const;

  const std::vector<T>& Derivatives() const {
    return mDerivatives;
  }

protected:
  std::pair<T, T> Interpolate(const double time, const int lowerIdx) const;

private:
  // for efficiency
  CubicHermiteInterpolator(
      const std::pair<std::vector<double>, std::vector<T>>& vectors);

  std::vector<T> mDerivatives;
};

#endif // SKYNET_UTILITIES_CUBICHERMITEINTERPOLATOR_HPP_
