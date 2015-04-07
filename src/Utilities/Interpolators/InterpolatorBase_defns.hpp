/// \file InterpolatorBase.cpp
/// \author jlippuner
/// \since Jan 7, 2015
///
/// \brief
///
///

#include "Utilities/Interpolators/InterpolatorBase.hpp"

#include <algorithm>
#include <stdexcept>

namespace { // unnamed so this function can only be used in this file

struct CompareTimes: std::binary_function<std::size_t, std::size_t, bool> {
  CompareTimes(const std::vector<double>& times) :
      mTimes(times) {}

  bool operator()(std::size_t lhs, std::size_t rhs) const {
    return mTimes[lhs] < mTimes[rhs];
  }

  const std::vector<double>& mTimes;
};

} // namespace [unnamed]

template<typename T>
InterpolatorBase<T>::InterpolatorBase(const std::vector<double>& times,
    const std::vector<T>& values) {
  if (times.size() != values.size())
    throw std::invalid_argument("Attempted to create an InterpolatorBase "
        "with different numbers of times and values");

  if (times.size() < 2)
    throw std::invalid_argument("Attempted to create an InterpolatorBase "
        "with less than 2 times");

  // get time sorting order
  std::vector<std::size_t> order(times.size());
  for (unsigned int i = 0; i < times.size(); ++i)
    order[i] = i;

  std::stable_sort(order.begin(), order.end(), CompareTimes(times));

  mTimes.reserve(times.size());
  mValues.reserve(times.size());

  mTimes.push_back(times[order[0]]);
  mValues.push_back(values[order[0]]);

  for (unsigned int i = 1; i < order.size(); ++i) {
    auto lastIdx = mTimes.size() - 1;

    if (times[order[i]] > mTimes[lastIdx]) {
      mTimes.push_back(times[order[i]]);
      mValues.push_back(values[order[i]]);
    } else {
      mValues[lastIdx] = values[order[i]];
    }
  }
}

template<typename T>
std::pair<T, T> InterpolatorBase<T>::ValueAndFirstDeriv(
    const double time) const {
  if (time < MinTime())
    throw std::invalid_argument("Tried to interpolate to a time that is "
        "smaller than the minimum known time (" + std::to_string(time) + " < "
        + std::to_string(MinTime()) + ").");
  if (time > MaxTime())
    throw std::invalid_argument("Tried to interpolate to a time that is larger "
        "than the maximum known time (" + std::to_string(time) + " > "
        + std::to_string(MaxTime()) + ").");

  if (time == mTimes[0])
    return Interpolate(time, 0);

  for (unsigned int i = 1; i < mTimes.size() - 1; ++i) {
    if (time == mTimes[i])
      return Interpolate(time, i);

    // TODO make this a bisection
    if (mTimes[i] > time)
      return Interpolate(time, i - 1);
  }

  return Interpolate(time, mTimes.size() - 2);
}

template<typename T>
T InterpolatorBase<T>::operator()(const double time) const {
  return ValueAndFirstDeriv(time).first;
}

template<typename T>
T InterpolatorBase<T>::FirstDerivative(const double time) const {
  return ValueAndFirstDeriv(time).second;
}

template<typename T>
std::pair<std::vector<double>, std::vector<T>>
InterpolatorBase<T>::SplitUpPoints(
    const std::vector<std::pair<double, T>>& points) {
  std::vector<double> times(points.size());
  std::vector<T> values(points.size());

  for (unsigned int i = 0; i < points.size(); ++i) {
    times[i] = points[i].first;
    values[i] = points[i].second;
  }

  return {times, values};
}
