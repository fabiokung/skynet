/// \file PiecewiseLinearFunction.cpp
/// \author jlippuner
/// \since Jul 8, 2014
///
/// \brief
///
///

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <valarray>

#include "Utilities/CodeError.hpp"

namespace { // unnamed so this function can only be used in this file

template<typename T>
bool HasZeroValues(const T& value);

template<>
bool HasZeroValues<double>(const double& value) {
  return value <= 0;
}

template<>
bool HasZeroValues<std::valarray<double>>(const std::valarray<double>& value) {
  auto isZero = (value <= 0.0);
  return isZero.max();
}

} // namespace [unnamed]

template<typename T>
GeneralPiecewiseLinearFunction<T>::GeneralPiecewiseLinearFunction(
    const std::vector<double>& times, const std::vector<T>& values,
    const bool interpolateInLogSpace) :
    InterpolatorBase<T>(times, values),
    mInterpolationInLogSpace(interpolateInLogSpace),
    mOffset(0.0) {
  // we can only assign the adjusted values and times now, because the
  // constructor of the InterpolatorBase might change the order
  mAdjustedTimes = this->Times();
  mAdjustedValues = this->Values();

  if (mInterpolationInLogSpace) {
    if (this->MinTime() < 1.0)
      mOffset = 1.0 - this->MinTime();

    for (unsigned int i = 0; i < this->Times().size(); ++i) {
      mAdjustedTimes[i] = log(this->Times()[i] + mOffset);

      //if (mPoints[i].second <= 0.0)
      if (HasZeroValues<T>(this->Values()[i]))
        throw std::invalid_argument("Tried to do interpolation in log space "
            "with negative y values");

      mAdjustedValues[i] = log(this->Values()[i]);
    }
  }
}

template<typename T>
GeneralPiecewiseLinearFunction<T>::GeneralPiecewiseLinearFunction(
    const std::vector<std::pair<double, T>>& points,
    const bool interpolateInLogSpace) :
    GeneralPiecewiseLinearFunction<T>(
        InterpolatorBase<T>::SplitUpPoints(points), interpolateInLogSpace) {}

template<typename T>
GeneralPiecewiseLinearFunction<T>::GeneralPiecewiseLinearFunction(
    const std::pair<std::vector<double>, std::vector<T>>& vectors,
    const bool interpolateInLogSpace) :
    GeneralPiecewiseLinearFunction<T>(vectors.first, vectors.second,
        interpolateInLogSpace) {}

template<typename T>
std::unique_ptr<FunctionVsTime<T>>
GeneralPiecewiseLinearFunction<T>::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<T>>(
      new GeneralPiecewiseLinearFunction<T>(*this));
}

template<typename T>
std::pair<T, T> GeneralPiecewiseLinearFunction<T>::Interpolate(
    const double time, const int lowerIdx) const {
  double timeLow = mAdjustedTimes[lowerIdx];
  double timeHigh = mAdjustedTimes[lowerIdx + 1];
  T valueLow = mAdjustedValues[lowerIdx];
  T valueHigh = mAdjustedValues[lowerIdx + 1];

  double t = time;
  if (mInterpolationInLogSpace)
    t = log(time + mOffset);

  T firstDeriv = (valueHigh - valueLow) / (timeHigh - timeLow);
  T value = valueLow + firstDeriv * (t - timeLow);

  if (mInterpolationInLogSpace) {
    value = exp(value);
    firstDeriv *= value / (time + mOffset);
  }

  return {value, firstDeriv};
}

// explicit template instantiations
template class GeneralPiecewiseLinearFunction<double>;
template class GeneralPiecewiseLinearFunction<std::valarray<double>>;
