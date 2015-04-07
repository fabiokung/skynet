/// \file GeneralRungeKutta.hpp
/// \author jlippuner
/// \since Sep 16, 2014
///
/// \brief General Runge-Kutta integrator.
///
///

#ifndef SKYNET_GENERALRUNGEKUTTA_HPP_
#define SKYNET_GENERALRUNGEKUTTA_HPP_

#include <array>
#include <cassert>
#include <functional>
#include <vector>

#include "Utilities/FloatingPointComparison.hpp"

/// This is a general Runge-Kutta integrator of order ORDER. It operates on the
/// data type T, which must have the operators +, +=, and * overloaded. The
/// right-hand side function takes as arguments the current time and the current
/// value of the integrated variable. The Runge-Kutta method is determined by
/// specifying the a_{ij}, b_i, and c_i coefficients of the Butcher tableau (see
/// http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge.E2.80.93Kutta_methods).
/// NUM sets of b_i coefficients are required, so that NUM different results
/// can be reported. This is useful for combining different orders of RK methods
/// that share the same a_{ij} and c_i coefficients (e.g. Runge-Kutta-Fehlberg
/// method).
template<typename T, int ORDER, int NUM>
class GeneralRungeKutta {
  static_assert(ORDER >= 1, "The order of a Runge-Kutta method must be at "
      "lest 1.");
  static_assert(NUM >= 1, "A Runge-Kutta method must return at least 1 "
      "result.");
  public:
  GeneralRungeKutta(const std::function<T(const double, const T&)>& rhs,
      const std::array<std::vector<double>, ORDER - 1>& a,
      const std::array<std::array<double, ORDER>, NUM>& b,
      const std::array<double, ORDER - 1>& c) :
      mRhs(rhs),
      mA(a),
      mB(b),
      mC(c) {}

  /// Returns the NUM new values at time t + dt given the previous value prev at
  /// time t.
  std::array<T, NUM> NextStep(const T& prev, const double t,
      const double dt) const;

private:
  const std::function<T(const double, const T&)> mRhs;
  const std::array<std::vector<double>, ORDER - 1> mA;
  const std::array<std::array<double, ORDER>, NUM> mB;
  const std::array<double, ORDER - 1> mC;
};

template<typename T, int ORDER, int NUM>
std::array<T, NUM> GeneralRungeKutta<T, ORDER, NUM>::NextStep(const T& prev,
    const double t, const double dt) const {
  std::array<T, ORDER> ks;

  // first k is special
  ks[0] = mRhs(t, prev) * dt;

  for (int i = 1; i < ORDER; ++i) {
    assert((int)mA[i - 1].size() == i);

    // for consistency check
    double sumA = 0.0;

    T argument = prev;

    for (int j = 0; j < i; ++j) {
      sumA += mA[i - 1][j];
      argument += ks[j] * mA[i - 1][j];
    }

    // check consistency
    assert(fcmp(sumA, mC[i - 1]) == 0);

    ks[i] = mRhs(t + mC[i - 1] * dt, argument) * dt;
  }

  std::array<T, NUM> res;

  for (int n = 0; n < NUM; ++n) {
    res[n] = prev;
    for (int i = 0; i < ORDER; ++i)
      res[n] += ks[i] * mB[n][i];
  }

  return res;
}

#endif // SKYNET_GENERALRUNGEKUTTA_HPP_
