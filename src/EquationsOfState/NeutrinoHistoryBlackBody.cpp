/// \file NeutrinoHistoryBlackBody.cpp
/// \author jlippuner
/// \since Dec 22, 2015
///
/// \brief
///
///

#include <gsl/gsl_sf_fermi_dirac.h>

#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"

std::shared_ptr<NeutrinoDistribution> NeutrinoHistoryBlackBody::operator()(
    const double time) const {
  auto T9 = mT9VsTime(time);
  auto eta = mEtaVsTime(time);
  auto Lnu = mLVsTime(time);

  // need cast to double to avoid linking error in debug mode
  std::valarray<double> TInMeV = T9
      * (double)Constants::BoltzmannConstantInMeVPerGK;

  std::valarray<double> norm(mSpecies.size());
  std::valarray<double> BBconst(mSpecies.size());
  std::valarray<double> rnu(mSpecies.size());
  double radius = mRadiusVsTime(time);

  for (unsigned int i = 0; i < mSpecies.size(); i++) {
    // The energy flux of a pinched Fermi-Dirac spectrum scales with F_3(eta),
    // so the blackbody constant must be eta-dependent; using the eta = 0 value
    // (7 pi^3 / 240) overestimates the number flux for eta > 0 (6.6x at
    // eta = 2.1, 3.9x at eta = 1.5).
    BBconst[i] = 0.5 / Constants::Pi * tgamma(4.0)
        * gsl_sf_fermi_dirac_int(3, eta[i])
        * Constants::ErgPerMeV
        / pow(Constants::ReducedPlanckConstantInMeVSec, 3)
        / pow(Constants::SpeedOfLightInCmPerSec, 2);
    rnu[i] = sqrt(Lnu[i] / BBconst[i] / pow(TInMeV[i], 4.0));

    // This function returns the integration over solid angle for a distribution
    // function with a minimum mu cutoff normalized by four pi.  This is the
    // expectation for blackbody emitting sphere.
    if (mPointSource) {
      norm[i] = Lnu[i] / (pow(TInMeV[i], 4) * pow(radius, 2) * 4.0 * BBconst[i]);
    } else if (rnu[i] * rnu[i] <= radius * radius) {
      norm[i] = 0.5 - 0.5 * sqrt(1.0 - rnu[i] * rnu[i] / (radius * radius));
    } else {
    // Model the neutrino decoupling region with an exponential decrease.
    // This model makes little physical sense but should smoothly
    // transition from equilibrium to the free streaming limit
      norm[i] = 1.0
          - 0.5 * exp(-6.0 * (rnu[i] * rnu[i] / (radius * radius) - 1.0));
    }
  }

  return NeutrinoDistributionFermiDirac::Create(mSpecies, T9, eta, norm);
}

void NeutrinoHistoryBlackBody::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# Neutrino History: NeutrinoHistoryThermal\n");
  pOutput->Log("#   Type: %s\n", mIsConst ? "constant" : "time-dependent");
  pOutput->Log("#   initial radius: %.5e\n", mRadiusVsTime.Values()[0]);
  pOutput->Log("#   Number of distributions: %lu\n", mSpecies.size());

  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    pOutput->Log("#   Distribution %u:\n", i);
    pOutput->Log("#     Species: %s\n",
        NeutrinoDistribution::NeutrinoSpeciesToStr(mSpecies[i]).c_str());
    pOutput->Log("#     initial temperature: %.5e\n", mT9VsTime.Values()[0][i]);
    pOutput->Log("#     initial luminosity: %.5e\n", mLVsTime.Values()[0][i]);
  }
}


