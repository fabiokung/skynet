/// \file NeutrinoDistributionTabulated.cpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#include "EquationsOfState/NeutrinoDistributionTabulated.hpp"

#include <stdexcept>

std::shared_ptr<NeutrinoDistribution> NeutrinoDistributionTabulated::Create(
    const std::valarray<NeutrinoSpecies>& species,
    const std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>>&
        spectra,
    const std::valarray<double>& normalization) {
  return std::shared_ptr<NeutrinoDistribution>(
      new NeutrinoDistributionTabulated(species, spectra, normalization));
}

void NeutrinoDistributionTabulated::SetLocalT9(const double T9) {
  NeutrinoDistribution::SetLocalT9(T9);

  for (const auto& spectrum : mSpectra) {
    if (mLocalTMeV > spectrum->TailTemperatureMeV())
      throw std::runtime_error("Local fluid temperature "
          + std::to_string(mLocalTMeV) + " MeV exceeds the tabulated neutrino "
          "spectrum's effective temperature "
          + std::to_string(spectrum->TailTemperatureMeV()) + " MeV. The "
          "Fermi-Dirac path would thermalize the spectrum up to the local "
          "temperature here (the trapping floor); a tabulated spectrum cannot, "
          "so the two diverge. Use this spectrum only where the local fluid "
          "stays cooler than it.");
  }
}

std::function<double(double)>
NeutrinoDistributionTabulated::DistributionFunction(
    const NeutrinoSpecies species) const {
  std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>> matching;
  std::vector<double> norms;

  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    if (mSpecies[i] != species)
      continue;

    matching.push_back(mSpectra[i]);
    norms.push_back(mNormalization[i]);
  }

  if (matching.empty())
    return [] (const double /*enuInMeV*/) { return 0.0; };

  if (matching.size() == 1) {
    auto spectrum = matching[0];
    const double norm = norms[0];

    return [=] (const double enuInMeV) {
      return norm * (*spectrum)(enuInMeV);
    };
  }

  return [=] (const double enuInMeV) {
    double sum = 0.0;
    for (unsigned int i = 0; i < matching.size(); ++i)
      sum += norms[i] * (*matching[i])(enuInMeV);
    return sum;
  };
}

NeutrinoDistributionTabulated::NeutrinoDistributionTabulated(
    const std::valarray<NeutrinoSpecies>& species,
    const std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>>&
        spectra,
    const std::valarray<double>& normalization) :
    NeutrinoDistribution(species),
    mSpectra(spectra),
    mNormalization(normalization) {
  if (species.size() != spectra.size())
    throw std::invalid_argument("Tabulated neutrino distribution has "
        + std::to_string(species.size()) + " species but "
        + std::to_string(spectra.size()) + " spectra");

  if (species.size() != normalization.size())
    throw std::invalid_argument("Tabulated neutrino distribution has "
        + std::to_string(species.size()) + " species but "
        + std::to_string(normalization.size()) + " normalizations");
}
