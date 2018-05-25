/// \file Neutrino.cpp
/// \author lroberts
/// \since May 13, 2015
///
/// \brief
///
///

#include "Reactions/Neutrino.hpp"

#include "Reactions/REACLIB.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Network/NetworkOptions.hpp"

Neutrino::Neutrino(const std::string& pathToFile,
    const NuclideLibrary& nuclib) :
    mSource(pathToFile) {
  // use the same file format as REACLIB
  REACLIB inp(pathToFile);

  auto entries = inp.Entries();

  for (auto& entry : entries) {
    if (REACLIB::ReactionContainsUnkownNuclides(entry, nuclib))
      continue;

    // the first rate fitting coefficient is the matrix element and the second
    // is the weak magnetism correction, all other coefficients must be 0
    if (entry.RateFittingCoefficients().size() < 2) {
      throw std::runtime_error("Got an invalid neutrino reaction. It must "
          "have at least 2 rate fitting coefficients.");
    }

    for (unsigned i = 2; i < entry.RateFittingCoefficients().size(); ++i) {
      if (entry.RateFittingCoefficients()[i] != 0.0) {
        throw std::runtime_error("Got an invalid neutrino reaction. It has "
            "more than 2 non-zero rate fitting coefficients");
      }
    }

    std::vector<int> reactantAs, productAs, reactantZs, productZs;
    int totalReactantZ = 0;
    int totalProductZ = 0;

    for (auto react : entry.Reactants()) {
      int id = nuclib.NuclideIdsVsNames().at(react);
      reactantAs.push_back(nuclib.As()[id]);
      reactantZs.push_back(nuclib.Zs()[id]);

      totalReactantZ += nuclib.Zs()[id];
    }

    for (auto react : entry.Products()) {
      int id = nuclib.NuclideIdsVsNames().at(react);
      productAs.push_back(nuclib.As()[id]);
      productZs.push_back(nuclib.Zs()[id]);

      totalProductZ += nuclib.Zs()[id];
    }

    if (abs(totalReactantZ - totalProductZ) != 1) {
      throw std::runtime_error("Got an invalid neutrino reaction. Total Z of "
          "reactants and products must differ by exactly 1");
    }

    double matrixElement = entry.RateFittingCoefficients()[0];
    double Wm = entry.RateFittingCoefficients()[1];

    mEntries.push_back(NeutrinoEntry(reactantAs, reactantZs, productAs,
        productZs, entry.Q(), matrixElement, Wm,
        (totalReactantZ >= totalProductZ)));
  }
}

Neutrino Neutrino::FromREACLIBDecays(const std::string& REACLIBfile,
    const NuclideLibrary& nuclib, NetworkOptions opts) {

  REACLIBReactionLibrary weakReactionLibrary(REACLIBfile,
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "REACLIB Weak reactions", nuclib, opts, false);

  std::vector<NeutrinoEntry> nuRates;

  for (unsigned int i=0; i<weakReactionLibrary.Reactions().size(); ++i) {
    std::cout << weakReactionLibrary.Reactions()[i].String() << " "
        << exp(weakReactionLibrary.RateFittingCoefficients()[0][i]) << std::endl;
    if (!weakReactionLibrary.RateIsTemperatureDependent()[i]) {
      try { // This try statement should prevent any reactions containing
            // more than one product and reactant from being added, since
            // FromVacuumBetaDecay throws an error in this case.
      // Forward rate
      nuRates.push_back(NeutrinoEntry::FromVacuumBetaDecay(
          weakReactionLibrary.Reactions()[i],
          exp(weakReactionLibrary.RateFittingCoefficients()[0][i]),
          nuclib, false));
      // Reverse rate corresponding to beta decay
      nuRates.push_back(NeutrinoEntry::FromVacuumBetaDecay(
          weakReactionLibrary.Reactions()[i],
          exp(weakReactionLibrary.RateFittingCoefficients()[0][i]),
          nuclib, true));
      } catch(...) {}
    }
  }

  return Neutrino(nuRates, REACLIBfile);
}

std::vector<Reaction> Neutrino::GetValidReactions(
    const NuclideLibrary& nuclib) const {
  std::vector<Reaction> reactions;
  for (auto entry : mEntries) {
    auto reac = entry.GetReaction(nuclib);
    if (reac.ContainsOnlyKnownNuclides(nuclib))
      reactions.push_back(reac);
  }
  return reactions;
}
