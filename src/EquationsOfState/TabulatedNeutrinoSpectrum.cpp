/// \file TabulatedNeutrinoSpectrum.cpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#include "EquationsOfState/TabulatedNeutrinoSpectrum.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace { // unnamed so these can only be used in this file

// 8-point Gauss-Legendre on [-1, 1]
const double GLNodes[4] = { 0.1834346424956498, 0.5255324099163290,
    0.7966664774136267, 0.9602898564975363 };
const double GLWeights[4] = { 0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763 };

// exp() traps on overflow (FE_OVERFLOW is enabled). The constructor rejects
// any knot or below-grid value with log f > MaxLog, so the cap in operator()
// only ever guards a cubic-Hermite overshoot between knots.
const double MaxLog = 700.0;

// largest share of the energy flux allowed above the top of the energy grid
const double MaxTailFraction = 1.0e-2;

void ValidateGrid(const std::vector<double>& energiesMeV,
    const std::size_t numValues) {
  if (energiesMeV.size() != numValues)
    throw std::invalid_argument("Tabulated neutrino spectrum has "
        + std::to_string(energiesMeV.size()) + " energies but "
        + std::to_string(numValues) + " values");

  if (energiesMeV.size() < 3)
    throw std::invalid_argument("Tabulated neutrino spectrum needs at least "
        "3 energy points");

  for (unsigned int i = 0; i < energiesMeV.size(); ++i) {
    if (i > 0 && energiesMeV[i] <= energiesMeV[i-1])
      throw std::invalid_argument("Tabulated neutrino spectrum energies must "
          "be strictly increasing, but E[" + std::to_string(i) + "] = "
          + std::to_string(energiesMeV[i]) + " <= E[" + std::to_string(i-1)
          + "] = " + std::to_string(energiesMeV[i-1]));

    if (energiesMeV[i] < 0.0)
      throw std::invalid_argument("Tabulated neutrino spectrum has a negative "
          "energy " + std::to_string(energiesMeV[i]) + " MeV");
  }
}

std::vector<double> TakeLogs(const std::vector<double>& energiesMeV,
    const std::vector<double>& values) {
  ValidateGrid(energiesMeV, values.size());

  std::vector<double> logs(values.size());

  for (unsigned int i = 0; i < values.size(); ++i) {
    if (!(values[i] > 0.0) || !std::isfinite(values[i]))
      throw std::invalid_argument("Tabulated neutrino spectrum must be "
          "positive and finite, but f(" + std::to_string(energiesMeV[i])
          + " MeV) = " + std::to_string(values[i]));

    logs[i] = log(values[i]);

    if (logs[i] > MaxLog)
      throw std::invalid_argument("Tabulated neutrino spectrum value f("
          + std::to_string(energiesMeV[i]) + " MeV) = "
          + std::to_string(values[i]) + " is too large to exponentiate "
          "(log f = " + std::to_string(logs[i]) + " > " + std::to_string(MaxLog)
          + "); occupation numbers are of order 1");
  }

  return logs;
}

// int_0^infinity (u + e0)^n exp(-b u) du, the E > MaxEnergy tail
double TailIntegral(const int n, const double e0, const double b) {
  double binom = 1.0;
  double factorial = 1.0;
  double bPow = b;
  double sum = 0.0;

  for (int k = 0; k <= n; ++k) {
    if (k > 0) {
      binom = binom * (n - k + 1) / k;
      factorial *= k;
      bPow *= b;
    }
    sum += binom * pow(e0, n - k) * factorial / bPow;
  }

  return sum;
}

} // namespace [unnamed]

TabulatedNeutrinoSpectrum::TabulatedNeutrinoSpectrum(
    const std::vector<double>& energiesMeV,
    const std::vector<double>& values) :
    mLogF(energiesMeV, TakeLogs(energiesMeV, values)) {
  mLogFLow = mLogF.Values().front();
  mLogFHigh = mLogF.Values().back();
  mSlopeLow = mLogF.FirstDerivative(mLogF.MinTime());
  mSlopeHigh = mLogF.FirstDerivative(mLogF.MaxTime());

  if (!(mSlopeHigh < 0.0))
    throw std::invalid_argument("Tabulated neutrino spectrum does not decay at "
        "the top of the energy grid (d log f / dE = "
        + std::to_string(mSlopeHigh) + " at " + std::to_string(mLogF.MaxTime())
        + " MeV); the E -> infinity tail is not integrable");

  // The rate integrals query down to E = 0; below the grid f is continued
  // log-linearly, which for a decaying spectrum grows toward low E. Reject a
  // grid whose extrapolation to E = 0 would overflow exp().
  const double logFAtZero = mLogFLow + mSlopeLow * (0.0 - MinEnergy());
  if (logFAtZero > MaxLog)
    throw std::invalid_argument("Tabulated neutrino spectrum extrapolates below "
        "its grid (starts at " + std::to_string(MinEnergy()) + " MeV) to "
        "log f(0) = " + std::to_string(logFAtZero) + " > " + std::to_string(MaxLog)
        + ", too large to exponentiate; start the grid nearer E = 0");

  // Everything above the grid is a straight line in log f, so the grid has to
  // reach far enough into the tail that the extrapolation carries a negligible
  // part of the energy flux. Past a percent the spectrum is set by the
  // extrapolation rather than by the table.
  const double tailFraction = TailMoment(3) / EnergyMoment(3);
  if (tailFraction > MaxTailFraction)
    throw std::invalid_argument("Tabulated neutrino spectrum carries "
        + std::to_string(100.0 * tailFraction) + "% of its energy flux above "
        "the top of the energy grid (" + std::to_string(mLogF.MaxTime())
        + " MeV), where it is extrapolated; extend the grid");
}

double TabulatedNeutrinoSpectrum::operator()(const double eNuMeV) const {
  double logF;

  if (eNuMeV < MinEnergy())
    logF = mLogFLow + mSlopeLow * (eNuMeV - MinEnergy());
  else if (eNuMeV > MaxEnergy())
    logF = mLogFHigh + mSlopeHigh * (eNuMeV - MaxEnergy());
  else
    logF = mLogF(eNuMeV);

  // the min is a backstop for inter-knot overshoot; the constructor has already
  // bounded the knots and the below-grid extrapolation to log f <= MaxLog
  return exp(std::min(logF, MaxLog));
}

double TabulatedNeutrinoSpectrum::IntervalMoment(const int n,
    const double eLow, const double eHigh) const {
  const double half = 0.5 * (eHigh - eLow);
  const double mid = 0.5 * (eHigh + eLow);
  double sum = 0.0;

  for (int i = 0; i < 4; ++i) {
    const double eMinus = mid - half * GLNodes[i];
    const double ePlus = mid + half * GLNodes[i];
    sum += GLWeights[i] * (pow(eMinus, n) * (*this)(eMinus)
        + pow(ePlus, n) * (*this)(ePlus));
  }

  return half * sum;
}

double TabulatedNeutrinoSpectrum::TailMoment(const int n) const {
  return exp(std::min(mLogFHigh, MaxLog))
      * TailIntegral(n, MaxEnergy(), -mSlopeHigh);
}

double TabulatedNeutrinoSpectrum::EnergyMoment(const int n) const {
  const auto& energies = mLogF.Times();
  double moment = 0.0;

  if (MinEnergy() > 0.0)
    moment += IntervalMoment(n, 0.0, MinEnergy());

  for (unsigned int i = 0; i + 1 < energies.size(); ++i)
    moment += IntervalMoment(n, energies[i], energies[i+1]);

  return moment + TailMoment(n);
}
