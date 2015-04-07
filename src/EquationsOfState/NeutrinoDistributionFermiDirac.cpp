/// \file NeutrinoDistributionFermiDirac.cpp
/// \author jlippuner
/// \since Feb 4, 2016
///
/// \brief
///
///

#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"

#include <cassert>

std::shared_ptr<NeutrinoDistribution> NeutrinoDistributionFermiDirac::Create(
    const std::valarray<NeutrinoSpecies>& species,
    const std::valarray<double>& T9, const std::valarray<double>& eta) {
  return std::shared_ptr<NeutrinoDistribution>(
      new NeutrinoDistributionFermiDirac(species, T9, eta,
          std::valarray<double>(1.0, species.size())));
}

std::shared_ptr<NeutrinoDistribution> NeutrinoDistributionFermiDirac::Create(
    const std::valarray<NeutrinoSpecies>& species,
    const std::valarray<double>& T9, const std::valarray<double>& eta,
    const std::valarray<double>& normalization) {
  return std::shared_ptr<NeutrinoDistribution>(
      new NeutrinoDistributionFermiDirac(species, T9, eta, normalization));
}

std::function<double(double)>
NeutrinoDistributionFermiDirac::DistributionFunction(
    const NeutrinoSpecies species) const {
  // count how many species match the given species
  int cnt = 0;
  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    if (mSpecies[i] == species)
      ++cnt;
  }

  // there are no matching neutrino species, hence the distribution function
  // is 0.0
  if (cnt == 0)
    return [] (const double /*enuInMeV*/) { return 0.0; };

  // now collect temperatures and chemical potential of the matching species
  std::valarray<double> Ts(cnt);
  std::valarray<double> etas(cnt);
  std::valarray<double> norms(cnt);

  int idx = 0;
  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    if (mSpecies[i] != species)
      continue;

    if (mTInMeV[i] < mLocalTMeV)
      Ts[idx] = mLocalTMeV;
    else
      Ts[idx] = mTInMeV[i];

    etas[idx] = mEta[i];

    norms[idx] = mNormalization[i];

    ++idx;
  }

  if (cnt == 1) {
    double T = Ts[0];
    double eta = etas[0];
    double norm = norms[0];

    return [=] (const double enuInMeV) {
      double arg = enuInMeV / T - eta;
      if (arg > 100.0)
        arg = 100.0;
      return norm / (exp(arg) + 1.0);
    };
  } else {
    return [=] (const double enuInMeV) {
      std::valarray<double> enu(enuInMeV, Ts.size());
      std::valarray<double> arg = enu / Ts - etas;

      for (unsigned int i = 0; i < arg.size(); ++i)
        arg[i] = std::min(arg[i], 100.0);

      std::valarray<double> fd = norms / (exp(arg) + 1.0);
      return fd.sum();
    };
  }
}

NeutrinoDistributionFermiDirac::NeutrinoDistributionFermiDirac(
    const std::valarray<NeutrinoSpecies>& species,
    const std::valarray<double>& T9, const std::valarray<double>& eta,
    const std::valarray<double>& normalization) :
    NeutrinoDistribution(species),
    // need cast to double to avoid linking error in debug mode
    mTInMeV((double)Constants::BoltzmannConstantInMeVPerGK * T9),
    mEta(eta),
    mNormalization(normalization) {
  assert(species.size() == T9.size());
  assert(species.size() == eta.size());
  assert(species.size() == normalization.size());
}
