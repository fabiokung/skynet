/// \file ConstReactionLibrary.cpp
/// \author jlippuner
/// \since Jan 13, 2015
///
/// \brief
///
///

#include "Reactions/ConstReactionLibrary.hpp"

ConstReactionLibrary::ConstReactionLibrary(const ReactionType type,
    const std::string& description, const std::string& source,
    const std::vector<Reaction>& reactions, const std::vector<double> rates,
    const NuclideLibrary& nuclib, const NetworkOptions& opts) :
    ReactionLibraryBase(type, description, source, reactions, nuclib, opts,
        false),
    mInverseRates() {
  std::vector<double> ratesWithFactorials(NumAllReactions());

  if (rates.size() != reactions.size())
    throw std::invalid_argument("Attempted to create a ConstReactionLibrary "
        "with different numbers of reactions and rates.");

  // divide rates by factorial product of reactants
  for (unsigned int i = 0; i < NumAllReactions(); ++i)
    ratesWithFactorials[i] = rates[i]
        / (double)Factorials()[i].NsFactorialProdcutOfReactants;

  mAllRates = ratesWithFactorials;

  PopulateRates();
}

void ConstReactionLibrary::DoLoopOverReactionData(
    const std::function<void(GeneralReactionData * const)>& func) {
  func(&mAllRates);

  PopulateRates();
}

void ConstReactionLibrary::PopulateRates() {
  mRates = mAllRates.ActiveData();
}
