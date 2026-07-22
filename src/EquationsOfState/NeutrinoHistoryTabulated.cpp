/// \file NeutrinoHistoryTabulated.cpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#include "EquationsOfState/NeutrinoHistoryTabulated.hpp"

std::shared_ptr<NeutrinoDistribution> NeutrinoHistoryTabulated::operator()(
    const double time) const {
  auto cached = std::atomic_load(&mCache);
  if (cached && cached->Time == time)
    return cached->Distribution;

  auto logValues = mLogSpectrumVsTime(time);
  auto Lnu = mLVsTime(time);
  const double radius = mRadiusVsTime(time);

  // BBconst with Gamma(4) * F_3(eta) * T^4 factored out, i.e. the conversion
  // from the third energy moment of the spectrum to an energy flux
  const double fluxConst = 0.5 / Constants::Pi * Constants::ErgPerMeV
      / pow(Constants::ReducedPlanckConstantInMeVSec, 3)
      / pow(Constants::SpeedOfLightInCmPerSec, 2);

  const std::size_t numEnergies = mEnergiesMeV.size();

  std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>> spectra(
      mSpecies.size());
  std::valarray<double> norm(mSpecies.size());

  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    std::vector<double> snapshot(numEnergies);
    for (unsigned int e = 0; e < numEnergies; ++e)
      snapshot[e] = exp(logValues[i * numEnergies + e]);

    std::shared_ptr<const TabulatedNeutrinoSpectrum> spectrum(
        new TabulatedNeutrinoSpectrum(mEnergiesMeV, snapshot));
    spectra[i] = spectrum;

    const double flux = fluxConst * spectrum->EnergyMoment(3);
    const double rnu2 = Lnu[i] / flux;

    // The point-source norm divides out EnergyMoment(3), so the spectrum's
    // overall scale cancels and only its shape matters. The sphere branches
    // below set the neutrinosphere radius from rnu2, which does depend on the
    // absolute third moment: there the spectrum must be a true occupation
    // number (max ~ 1), not just the right shape.
    if (mPointSource) {
      norm[i] = Lnu[i] / (4.0 * flux * radius * radius);
    } else if (rnu2 <= radius * radius) {
      norm[i] = 0.5 - 0.5 * sqrt(1.0 - rnu2 / (radius * radius));
    } else {
      // Model the neutrino decoupling region with an exponential decrease.
      // This model makes little physical sense but should smoothly
      // transition from equilibrium to the free streaming limit
      norm[i] = 1.0 - 0.5 * exp(-6.0 * (rnu2 / (radius * radius) - 1.0));
    }
  }

  auto distribution =
      NeutrinoDistributionTabulated::Create(mSpecies, spectra, norm);
  std::atomic_store(&mCache,
      std::shared_ptr<const Snapshot>(new Snapshot{ time, distribution }));

  return distribution;
}

void NeutrinoHistoryTabulated::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# Neutrino History: NeutrinoHistoryTabulated\n");
  pOutput->Log("#   Type: %s\n", mIsConst ? "constant" : "time-dependent");
  pOutput->Log("#   initial radius: %.5e\n", mRadiusVsTime.Values()[0]);
  pOutput->Log("#   Energy grid: %lu points, %.5e to %.5e MeV\n",
      mEnergiesMeV.size(), mEnergiesMeV.front(), mEnergiesMeV.back());
  pOutput->Log("#   Number of distributions: %lu\n", mSpecies.size());

  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    pOutput->Log("#   Distribution %u:\n", i);
    pOutput->Log("#     Species: %s\n",
        NeutrinoDistribution::NeutrinoSpeciesToStr(mSpecies[i]).c_str());
    pOutput->Log("#     initial luminosity: %.5e\n", mLVsTime.Values()[0][i]);
  }
}
