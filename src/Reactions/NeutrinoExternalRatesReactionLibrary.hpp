/// \file NeutrinoExternalRatesReactionLibrary.hpp
/// \author jlippuner
/// \since May 23, 2016
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_NEUTRINOEXTERNALRATESREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_NEUTRINOEXTERNALRATESREACTIONLIBRARY_HPP_

#include <vector>

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Neutrino.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

class NeutrinoExternalRatesReactionLibrary: public ReactionLibraryBase {
public:
  NeutrinoExternalRatesReactionLibrary(const Neutrino neutrinoLib,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, const std::vector<double>& times,
      const std::vector<double>& nuEAbsorptionRates,
      const std::vector<double>& antiNuEAbsorptionRates,
      const std::vector<double>& nuEEmissionRates,
      const std::vector<double>& antiNuEEmissionRates);

  NeutrinoExternalRatesReactionLibrary(const std::string& nuFile,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, const std::vector<double>& times,
      const std::vector<double>& nuEAbsorptionRates,
      const std::vector<double>& antiNuEAbsorptionRates,
      const std::vector<double>& nuEEmissionRates,
      const std::vector<double>& antiNuEEmissionRates) :
      NeutrinoExternalRatesReactionLibrary(Neutrino(nuFile, nuclib),
          description, nuclib, opts, times, nuEAbsorptionRates,
          antiNuEAbsorptionRates, nuEEmissionRates, antiNuEEmissionRates) {}

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new NeutrinoExternalRatesReactionLibrary(*this));
  }

  bool StateChanged(const ThermodynamicState /*thermoState*/,
      NetworkOutput * const /*pOutput*/) {
    // For now there is no state change, since the reactions
    // are assumed constant when they go off the edge of the
    // table.  At low temperature and density this is a good
    // assumption because capture rates go to zero and the
    // decay rates are the free space decay rates.  At high
    // density and temperature this is not so good, but if
    // we are out of the allowed upper limits of the tables
    // these reactions should definitely *not* be set to
    // zero.
    return false;
  }

  void SetTime(const double time) {
    mTime = time;
  }

  std::string Name() const {
    return "Neutrino Reaction External Rates Library";
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

  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

private:
  std::vector<double> mRates;
  std::vector<double> mInverseRates;

  // forward reaction is in the direction of electron or positron capture
  ReactionData<bool> mNueCap; // Is this an electron neutrino capture reaction

  double mTime;

  PiecewiseLinearFunction mNuEAbsRate, mAntiNuEAbsRate, mNuEEmRate,
      mAntiNuEEmRate;
};

#endif // SKYNET_REACTIONS_NEUTRINOEXTERNALRATESREACTIONLIBRARY_HPP_
