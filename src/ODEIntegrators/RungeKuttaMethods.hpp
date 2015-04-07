/// \file RungeKuttaMethods.hpp
/// \author jlippuner
/// \since Sep 16, 2014
///
/// \brief Some specific Runge-Kutte methods.
///
///

#ifndef SKYNET_RUNGEKUTTAMETHODS_HPP_
#define SKYNET_RUNGEKUTTAMETHODS_HPP_

#include <array>

#include "ODEIntegrators/GeneralRungeKutta.hpp"

/// The classical RK4.
template<typename T>
class RK4 {
public:
  explicit RK4(const std::function<T(const double, const T&)>& rhs) :
    mRK(rhs,
      {{
        {0.5},
        {0.0, 0.5},
        {0.0, 0.0, 1.0}
      }},
      {{
        {{1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0}}
      }},
      {{0.5, 0.5, 1.0}}) { }

  /// Returns the value at time t + dt given the previous value prev at time t.
  T NextStep(const T& prev, const double t, const double dt) const {
    return mRK.NextStep(prev, t, dt)[0];
  }

private:
  const GeneralRungeKutta<T, 4, 1> mRK;
};

/// The Runge-Kutta-Fehlberg method (RKF45).
template<typename T>
class RKF45 {
public:
  explicit RKF45(const std::function<T(const double, const T&)>& rhs) :
    mRK(rhs,
      {{
        {         0.25},
        {     3.0/32.0,       9.0/32.0},
        {1932.0/2197.0, -7200.0/2197.0,  7296.0/2197.0},
        {  439.0/216.0,           -8.0,   3680.0/513.0, -845.0/4104.0},
        {    -8.0/27.0,            2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0}
      }},
      {{
        {{16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0}},
        {{25.0/216.0, 0.0,  1408.0/2565.0,   2197.0/4104.0,  -1.0/5.0,      0.0}}
      }},
      {{0.25, 3.0/8.0, 12.0/13.0, 1.0, 0.5}}) { }

  /// Returns the value at time t + dt given the previous value prev at time t.
  std::array<T, 2> NextStep(const T& prev, const double t,
      const double dt) const {
    return mRK.NextStep(prev, t, dt);
  }

private:
  const GeneralRungeKutta<T, 6, 2> mRK;
};

#endif // SKYNET_RUNGEKUTTAMETHODS_HPP_
