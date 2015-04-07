/// \file ConstReactionLibrary.hpp
/// \author jlippuner
/// \since Jan 13, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_CONSTREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_CONSTREACTIONLIBRARY_HPP_

#include "Reactions/ReactionLibraryBase.hpp"

class ConstReactionLibrary: public ReactionLibraryBase {
public:
  ConstReactionLibrary(const ReactionType type, const std::string& description,
      const std::string& source, const std::vector<Reaction>& reactions,
      const std::vector<double> rates, const NuclideLibrary& nuclib,
      const NetworkOptions& opts);

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new ConstReactionLibrary(*this));
  }

  bool StateChanged(const ThermodynamicState /*thermoState*/,
      NetworkOutput * const /*pOutput*/) {
    // there is no state, so it never changes
    return false;
  }

  std::string Name() const {
    return "Constant Rates Reaction Library";
  }

  const std::vector<double>& Rates() const {
    return mRates;
  }

  const std::vector<double>& InverseRates() const {
    return mInverseRates;
  }

protected:
  void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func);

  void DoCalculateRates(const ThermodynamicState /*thermoState*/,
      const std::vector<double>& /*partitionFunctionsWithoutSpinTerms*/,
      const double /*expArgumentCap*/, const std::vector<int> * const /*pZs*/,
      const std::vector<double> * const /*pScreeningChemicalPotentialCorrection*/) {
    // do nothing
  }

private:
  void PopulateRates();

  ReactionData<double> mAllRates;

  std::vector<double> mRates;
  std::vector<double> mInverseRates;
};

#endif // SKYNET_REACTIONS_CONSTREACTIONLIBRARY_HPP_
