/// \file InterpolatorBase.hpp
/// \author jlippuner
/// \since Jan 7, 2015
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_INTERPOLATORS_INTERPOLATORBASE_HPP_
#define SKYNET_UTILITIES_INTERPOLATORS_INTERPOLATORBASE_HPP_

#include <vector>

#include "Utilities/FunctionVsTime.hpp"

template<typename T>
class InterpolatorBase : public FunctionVsTime<T> {
public:
  InterpolatorBase(const std::vector<double>& times,
      const std::vector<T>& values);

  // returns interpolated value and interpolated first derivative
  std::pair<T, T> ValueAndFirstDeriv(const double time) const;

  T operator()(const double time) const;

  T FirstDerivative(const double time) const;

  static std::pair<std::vector<double>, std::vector<T>> SplitUpPoints(
      const std::vector<std::pair<double, T>>& points);

  const std::vector<double>& Times() const {
    return mTimes;
  }

  const std::vector<T>& Values() const {
    return mValues;
  }

  double MinTime() const {
    return mTimes[0];
  }

  double MaxTime() const {
    return mTimes[mTimes.size()-1];
  }

protected:
  // This function must be implemented in the derived classes to actually do the
  // interpolation. When this function is called, it has already been ensured
  // that time is withing MinTime and MaxTime. lowerIdx is such that
  // mTimes[lowerIdx] <= time < mTimes[lowerIdx + 1], except if
  // lowerIdx == mTimes.size() - 2.
  // It returns the interpolated value and the first derivative.
  virtual std::pair<T, T> Interpolate(const double time,
      const int lowerIdx) const =0;

private:
  std::vector<double> mTimes;
  std::vector<T> mValues;
};

#include "Utilities/Interpolators/InterpolatorBase_defns.hpp"

#endif // SKYNET_UTILITIES_INTERPOLATORS_INTERPOLATORBASE_HPP_
