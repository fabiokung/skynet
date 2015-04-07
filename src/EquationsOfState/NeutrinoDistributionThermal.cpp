/// \file NeutrinoDistribution.cpp
/// \author lroberts
/// \since May 20, 2015
///
/// \brief
///
///

#include "EquationsOfState/NeutrinoDistributionThermal.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <algorithm>

#include "Utilities/Constants.hpp"

NeutrinoDistributionThermal::NeutrinoDistributionThermal(
    const std::valarray<NeutrinoSpecies>& species,
    const std::valarray<double>& T9, const std::valarray<double>& eta,
    const std::valarray<double>& deltaOmega, const double radCm,
    const std::valarray<double>& rnu) :
    mEta(eta),
    mTMeV(Constants::BoltzmannConstantInMeVPerGK * T9),
    mDeltaOmega(deltaOmega) {
  assert(species.size() == T9.size());
  assert(species.size() == eta.size());
  if (!deltaOmega.size() != 0)
    assert(species.size() == deltaOmega.size());
  if (!rnu.size() != 0)
    assert(species.size() == rnu.size());

  mSpecies = species;
  mRadius = radCm;
  mRnu = rnu;

  if (mDeltaOmega.size() == 0)
    mDeltaOmega = std::valarray<double>(1.0, species.size());
}

std::shared_ptr<NeutrinoDistribution>
NeutrinoDistributionThermal::FromTnuLnu(
    const std::valarray<NeutrinoSpecies>& species,
    const std::valarray<double>& T9, const std::valarray<double>& Lnu,
    const double radCm) {
  assert(species.size() == T9.size());
  assert(species.size() == Lnu.size());

  const double BBconst = 7.0 * pow(Constants::Pi, 3) / 240.0
      * Constants::ErgPerMeV
      / pow(Constants::ReducedPlanckConstantInMeVSec, 3)
      / pow(Constants::SpeedOfLightInCmPerSec, 2);

  std::valarray<double> eta(0.0, species.size());
  std::valarray<double> deltaOmega(species.size());

  std::valarray<double> TInMeV = T9 * Constants::BoltzmannConstantInMeVPerGK;
  std::valarray<double> rnu = sqrt(Lnu / BBconst / pow(TInMeV, 4.0));

  for (unsigned int i = 0; i < species.size(); i++) {
    deltaOmega[i] = NeutrinoDistributionThermal::GetAngularFactor(rnu[i],
        radCm);
  }

  return std::shared_ptr<NeutrinoDistribution>(
      new NeutrinoDistributionThermal(species, T9, eta, deltaOmega, radCm,
          rnu));
}



double NeutrinoDistributionThermal::GetAngularFactor(const double rnu,
    const double rad) {
  // This function returns the integration over solid angle for a distribution
  // function with a minimum mu cutoff normalized by four pi.  This is the
  // expectation for blackbody emitting sphere.
  if (rnu * rnu <= rad * rad) {
    return 0.5 - 0.5 * sqrt(1.0 - rnu * rnu / (rad * rad));
  }
  // Model the neutrino decoupling region with an exponential decrease.
  // This model makes little physical sense but should smoothly
  // transition from equilibrium to the free streaming limit
  return 1.0 - 0.5 * exp(-6.0 * (rnu * rnu / (rad * rad) - 1.0));
}

