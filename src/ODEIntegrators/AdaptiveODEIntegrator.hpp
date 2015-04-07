/// \file AdaptiveODEIntegrator.hpp
/// \author jlippuner
/// \since Sep 16, 2014
///
/// \brief
///
///

#ifndef SKYNET_ADAPTIVEODEITNEGRATOR_HPP_
#define SKYNET_ADAPTIVEODEITNEGRATOR_HPP_

#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <stdexcept>

/// This is the data type that the adaptive ODE integrator returns (and takes as
/// the argument for its TakeStep function). t is the current time, nextDt is
/// the time step to be used for the next step and value is the accepted value
/// of the current step.
template<typename T>
struct AdaptiveODEIntegratorResult {
  double t;
  double nextDt;
  T value;
};

/// This is a general adaptive ODE integrator that operates on a data type T. It
/// uses the function nextStep to call the ODE integrator to get the next step,
/// then uses the function error to calculate the error of this step. If the
/// error is larger than the maximum allowed error, the step is discarded, time
/// step is decreased and the step is attempted again with the reduced time
/// step. The time step is decreased (or increased) by a factor dtScaleFactor.
/// One also specifies a maximum and minimum allowed time step, as well as a
/// minimum error (when the error is below the minimum error, the time step
/// is increased).
template<typename T>
class AdaptiveODEIntegrator {
public:
  AdaptiveODEIntegrator(
      const std::function<T(const T&, const double, const double)>& nextStep,
      const std::function<double(const T&)>& error,
      const double dtMin, const double dtMax,
      const double errorMin, const double errorMax,
      const double dtScaleFactor) :
      mNextStep(nextStep),
      mError(error),
      mDtMin(dtMin), mDtMax(dtMax),
      mErrorMin(errorMin), mErrorMax(errorMax),
      mDtScaleFactor(dtScaleFactor) {
    assert(dtMin > 0.0);
    assert(dtMin < dtMax);
    assert(errorMin > 0.0);
    assert(errorMin < errorMax);
    assert(dtScaleFactor > 1.0);
  }

  /// Take a step and ensure that the error is within the tolerance.
  AdaptiveODEIntegratorResult<T> TakeStep(
      const AdaptiveODEIntegratorResult<T>& prev, const double dtMax) const;

private:
  const std::function<T(const T&, const double, const double)> mNextStep;
  const std::function<double(const T&)> mError;
  const double mDtMin, mDtMax;
  const double mErrorMin, mErrorMax;
  const double mDtScaleFactor;
};

template<typename T>
AdaptiveODEIntegratorResult<T> AdaptiveODEIntegrator<T>::TakeStep(
    const AdaptiveODEIntegratorResult<T>& prev, const double dtMax) const {
  double useDt = fmin(fmax(prev.nextDt, mDtMin), fmin(mDtMax, dtMax));

  // attempt step
  T next = mNextStep(prev.value, prev.t, useDt);
  double nextDt = useDt;
  double error = mError(next);

  if (error > mErrorMax) {
    // step is rejected, try again with smaller step size
    if (useDt <= mDtMin) {
      // the error is still too big, but we can't make the time step any
      // smaller, so just accept this step...

      fprintf(stderr, "WARNING: error = %e, dt = %e\n", error, useDt);
      nextDt = mDtMin;
    }
    else {
      AdaptiveODEIntegratorResult<T> newStep = prev;
      newStep.nextDt /= mDtScaleFactor;
      return TakeStep(newStep, dtMax);
    }
  }
  else {
    // we accept the step
    if (error < mErrorMin)
      // the error is below the tolerance, we can make the time step bigger
      nextDt = useDt * mDtScaleFactor;
    else
      // the error is within the tolerance, keep this time step
      nextDt = useDt;
  }

  AdaptiveODEIntegratorResult<T> res;
  res.t = prev.t + useDt;
  // make sure we are actually moving forward
  //assert (res.t > prev.t);

  res.nextDt = nextDt;
  res.value = next;

  return res;
}

#endif // SKYNET_ADAPTIVEODEITNEGRATOR_HPP_
