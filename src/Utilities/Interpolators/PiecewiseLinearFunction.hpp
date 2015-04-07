/// \file PiecewiseLinearFunction.hpp
/// \author jlippuner
/// \since Jul 8, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_PIECEWISELINEARFUNCTION_HPP_
#define SKYNET_UTILITIES_PIECEWISELINEARFUNCTION_HPP_

#include <vector>

#include "Utilities/Interpolators/InterpolatorBase.hpp"

template<typename T>
class GeneralPiecewiseLinearFunction : public InterpolatorBase<T> {
public:
  GeneralPiecewiseLinearFunction(const std::vector<double>& times,
      const std::vector<T>& values, const bool interpolateInLogSpace);

  GeneralPiecewiseLinearFunction(
      const std::vector<std::pair<double, T>>& points,
      const bool interpolateInLogSpace);

  std::unique_ptr<FunctionVsTime<T>> MakeUniquePtr() const;

protected:
  std::pair<T, T> Interpolate(const double time, const int lowerIdx) const;

private:
  // for efficiency
  GeneralPiecewiseLinearFunction(
      const std::pair<std::vector<double>, std::vector<T>>& vectors,
      const bool interpolateInLogSpace);

  // if the interpolation is not done in log space, these two vectors are
  // exactly the same as the original times and values, but if the interpolation
  // is done in log space, then these are the logs of the times (possibly
  // adjusted by the mOffset) and values
  std::vector<double> mAdjustedTimes;
  std::vector<T> mAdjustedValues;

  bool mInterpolationInLogSpace;

  // if the interpolation is done in log space, all x values must be positive,
  // to ensure this, the x values are offset by this amount internally
  double mOffset;
};

typedef GeneralPiecewiseLinearFunction<double> PiecewiseLinearFunction;

#endif // SKYNET_UTILITIES_PIECEWISELINEARFUNCTION_HPP_
