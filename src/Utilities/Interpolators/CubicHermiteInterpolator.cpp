/// \file CubicHermiteInterpolator.cpp
/// \author jlippuner
/// \since Jan 7, 2015
///
/// \brief
///
///

#include "Utilities/Interpolators/CubicHermiteInterpolator.hpp"

#include <valarray>

namespace { // unnamed so this function can only be used in this file

inline void TimeDiffs(const std::vector<double>& times, unsigned int i,
    double * const pD01, double * const pD02, double * const pD12) {
  double t0 = times[i-1];
  double t1 = times[i];
  double t2 = times[i+1];

  *pD01 = t0 - t1;
  *pD02 = t0 - t2;
  *pD12 = t1 - t2;
}

}

template<typename T>
CubicHermiteInterpolator<T>::CubicHermiteInterpolator(
    const std::vector<double>& times, const std::vector<T>& values) :
    InterpolatorBase<T>(times, values),
    mDerivatives(values.size()) {
  if (values.size() < 3)
    throw std::invalid_argument("Need at least 3 points for cubic Hermite "
        "interpolation");

  // calculate derivatives at all points, we use the 2nd order finite difference
  // formula obtained from the algorithm described in:
  // Bengt Fornberg, "Generation of Finite Difference Formulas on Arbitrarily
  // Spaced Grids" (1988), Mathematics of Computation 51:184 pp. 699-706
  //
  // for t[0] we use y[0], y[1], y[2]
  // for t[i] (0 < i < N-1) we use y[i-1], y[i], y[i+1]
  // for t[N-1] we use y[N-3], y[N-2], y[N-1]
  unsigned int N = values.size();

  // differences between x[i-1], x[i], and x[i+1]
  double d01, d02, d12;

  TimeDiffs(this->Times(), 1, &d01, &d02, &d12);

  // casting to T is necessary for the Intel compiler, otherwise the result will
  // be a 0-length valarray
  mDerivatives[0] = (T)(
      (d01 + d02) / (d01 * d02) * this->Values()[0]
          - d02 / (d01 * d12) * this->Values()[1]
          + d01 / (d02 * d12) * this->Values()[2]);

  for (unsigned int i = 1; i < N-1; ++i) {
    TimeDiffs(this->Times(), i, &d01, &d02, &d12);
    mDerivatives[i] = (T)(
        d12 / (d01 * d02) * this->Values()[i - 1]
            + (d01 - d12) / (d01 * d12) * this->Values()[i]
            - d01 / (d02 * d12) * this->Values()[i + 1]);
  }

  TimeDiffs(this->Times(), N-2, &d01, &d02, &d12);
  mDerivatives[N - 1] = (T)(
      -d12 / (d01 * d02) * this->Values()[N - 3]
          + d02 / (d01 * d12) * this->Values()[N - 2]
          - (d02 + d12) / (d02 * d12) * this->Values()[N - 1]);
}

template<typename T>
CubicHermiteInterpolator<T>::CubicHermiteInterpolator(
      const std::vector<std::pair<double, T>>& points) :
      CubicHermiteInterpolator<T>(InterpolatorBase<T>::SplitUpPoints(points)) {}

template<typename T>
CubicHermiteInterpolator<T>::CubicHermiteInterpolator(
    const std::pair<std::vector<double>, std::vector<T>>& vectors) :
    CubicHermiteInterpolator<T>(vectors.first, vectors.second) {}

template<typename T>
std::unique_ptr<FunctionVsTime<T>>
CubicHermiteInterpolator<T>::MakeUniquePtr() const {
  return std::unique_ptr<FunctionVsTime<T>>(
      new CubicHermiteInterpolator<T>(*this));
}

template<typename T>
std::pair<T, T> CubicHermiteInterpolator<T>::Interpolate(const double time,
    const int lowerIdx) const {
  double timeLow = this->Times()[lowerIdx];
  double timeHigh = this->Times()[lowerIdx + 1];
  T valueLow = this->Values()[lowerIdx];
  T valueHigh = this->Values()[lowerIdx + 1];
  T derivLow = mDerivatives[lowerIdx];
  T derivHigh = mDerivatives[lowerIdx + 1];

  double dt = timeHigh - timeLow;
  double t = (time - timeLow) / dt;
  double t2 = t * t;
  double t3 = t2 * t;

  T value = (2.0 * t3 - 3.0 * t2 + 1.0) * valueLow
      + (t3 - 2.0 * t2 + t) * dt * derivLow
      + (-2.0 * t3 + 3.0 * t2) * valueHigh
      + (t3 - t2) * dt * derivHigh;

  T firstDeriv = ((6.0 * t2 - 6.0 * t) * valueLow
      + (3.0 * t2 - 4.0 * t + 1.0) * dt * derivLow
      + (-6.0 * t2 + 6.0 * t) * valueHigh
      + (3.0 * t2 - 2.0 * t) * dt * derivHigh) / dt;

//  T secondDeriv = ((12.0 * t - 6.0) * valueLow
//      + (6.0 * t - 4.0) * dt * derivLow
//      + (-12.0 * t + 6.0) * valueHigh
//      + (6.0 * t - 2.0) * dt * derivHigh) / (dt * dt);
//
//  return {value, secondDeriv / (2.0 * firstDeriv)};
  return {value, firstDeriv};
}

// explicit template instantiations
template class CubicHermiteInterpolator<double>;
template class CubicHermiteInterpolator<std::valarray<double>>;
