/// \file TabulatedNeutrinoSpectrum.hpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_TABULATEDNEUTRINOSPECTRUM_HPP_
#define SRC_EQUATIONSOFSTATE_TABULATEDNEUTRINOSPECTRUM_HPP_

#include <vector>

#include "Utilities/Interpolators/CubicHermiteInterpolator.hpp"

// An unnormalized occupation number f(E) sampled on an arbitrary energy grid.
//
// The interpolation is a cubic Hermite in log f: exponentiating guarantees
// f > 0 at every query point, so a cubic can't ring the distribution negative
// the way it would in linear space (the network evaluates 1 - f(E) for Pauli
// blocking and integrates f over E to infinity, and neither survives a sign
// flip). For a pinched Fermi-Dirac spectrum log f is asymptotically linear in
// E, so the fit is near-exact in the tail: on a uniform 0.25 MeV grid to
// 100 MeV it reproduces T = 2.7 MeV, eta = 2.1 to 4e-6 pointwise and its third
// moment to 1e-10.
//
// Outside the grid f is continued log-linearly with the end-point slopes. The
// grid can be anything, but the spectrum has to be decaying at the top of it
// and has to have decayed far enough that the extrapolation carries no
// significant part of the E -> infinity tail the rate integrals run over.
class TabulatedNeutrinoSpectrum {
public:
  TabulatedNeutrinoSpectrum(const std::vector<double>& energiesMeV,
      const std::vector<double>& values);

  double operator()(const double eNuMeV) const;

  // int_0^infinity E^n f(E) dE, including the analytic log-linear tails
  double EnergyMoment(const int n) const;

  double MinEnergy() const {
    return mLogF.MinTime();
  }

  double MaxEnergy() const {
    return mLogF.MaxTime();
  }

  // Effective temperature in MeV from the high-energy tail: log f decays as
  // -E/T, so d(log f)/dE = -1/T. For a Fermi-Dirac shape this is exactly T.
  double TailTemperatureMeV() const {
    return -1.0 / mSlopeHigh;
  }

private:
  double IntervalMoment(const int n, const double eLow,
      const double eHigh) const;

  double TailMoment(const int n) const;

  CubicHermiteInterpolator<double> mLogF;
  double mLogFLow, mLogFHigh;   // log f at the grid ends
  double mSlopeLow, mSlopeHigh; // d(log f)/dE at the grid ends
};

#endif // SRC_EQUATIONSOFSTATE_TABULATEDNEUTRINOSPECTRUM_HPP_
