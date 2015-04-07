/// \file NeutrinoExternalRatesReactionLibrary
/// \author lroberts
/// \since May 23, 2016
///
/// \brief
///
///

#include "Reactions/NeutrinoExternalRatesReactionLibrary.hpp"
#include "Reactions/ReactionData.hpp"

#include <cmath>
#include <algorithm>

NeutrinoExternalRatesReactionLibrary::NeutrinoExternalRatesReactionLibrary(
    const Neutrino neutrinoLib, const std::string& description,
    const NuclideLibrary& nucLib, const NetworkOptions& opts,
    const std::vector<double>& times,
    const std::vector<double>& nuEAbsorptionRates,
    const std::vector<double>& antiNuEAbsorptionRates,
    const std::vector<double>& nuEEmissionRates,
    const std::vector<double>& antiNuEEmissionRates) :
    ReactionLibraryBase(ReactionType::Weak, description,
        neutrinoLib.GetSource(), neutrinoLib.GetValidReactions(nucLib),
        nucLib, opts, false),
    mTime(0.0),
    mNuEAbsRate(times, nuEAbsorptionRates, true),
    mAntiNuEAbsRate(times, antiNuEAbsorptionRates, true),
    mNuEEmRate(times, nuEEmissionRates, true),
    mAntiNuEEmRate(times, antiNuEEmissionRates, true) {

  auto entries = neutrinoLib.Entries();
  std::vector<bool> nueCap;
  for (auto entry : entries) {
    nueCap.push_back(entry.IsNueReaction());
  }

  mNueCap = ReactionData<bool>(nueCap);

  mRates = std::vector<double>(NumAllReactions());
  mInverseRates = std::vector<double>(NumAllReactions());
}

void NeutrinoExternalRatesReactionLibrary::DoLoopOverReactionData(
    const std::function<void(GeneralReactionData * const)>& func) {
  func(&mNueCap);
  mRates = std::vector<double>(Reactions().ActiveData().size());
  mInverseRates = std::vector<double>(Reactions().ActiveData().size());
}

void NeutrinoExternalRatesReactionLibrary::DoCalculateRates(
    const ThermodynamicState thermoState,
    const std::vector<double>& /*partitionFunctionsWithoutSpinTerms*/,
    const double /*expArgumentCap*/, const std::vector<int> * const /*pZs*/,
    const std::vector<double> * const /*pScreeningChemicalPotentialCorrection*/) {

  if (thermoState.T9() < Options().MinT9ForNeutrinoReactions) {
    std::fill(mRates.begin(), mRates.end(), 0.0);
    std::fill(mInverseRates.begin(), mInverseRates.end(), 0.0);
    return;
  }

  for (unsigned int i = 0; i < mNueCap.size(); i++) {
    NeutrinoSpecies nuSpec = NeutrinoSpecies::NuE;

    if (!mNueCap[i]) {
      nuSpec = NeutrinoSpecies::AntiNuE;
    }

    if (nuSpec == NeutrinoSpecies::NuE)
      mRates[i] = mNuEEmRate(mTime);
    else if (nuSpec == NeutrinoSpecies::AntiNuE)
      mRates[i] = mAntiNuEEmRate(mTime);
    else
      mRates[i] = 0.0;

    if (nuSpec == NeutrinoSpecies::NuE)
      mInverseRates[i] = mNuEAbsRate(mTime);
    else if (nuSpec == NeutrinoSpecies::AntiNuE)
      mInverseRates[i] = mAntiNuEAbsRate(mTime);
    else
      mInverseRates[i] = 0.0;
  }
}
