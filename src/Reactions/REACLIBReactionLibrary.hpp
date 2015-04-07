/// \file REACLIBReactionLibrary.hpp
/// \author jlippuner
/// \since Jul 24, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACLIBREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_REACLIBREACTIONLIBRARY_HPP_

#include <fstream>
#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/REACLIB.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

class REACLIBReactionLibrary: public ReactionLibraryBase {
public:
  void Dump(const std::string& fname) const {
    std::ofstream ofs(fname.c_str(),
        std::ostream::trunc | std::ostream::binary);
    boost::archive::binary_oarchive outArch(ofs);
    outArch << *this;
  }

  static REACLIBReactionLibrary ReadFromDisk(const std::string& fname,
      const NetworkOptions& opts) {
    REACLIBReactionLibrary lib(opts);
    std::ifstream ifs(fname.c_str(),
        std::ostream::in | std::ostream::binary);
    boost::archive::binary_iarchive inArch(ifs);
    inArch >> lib;
    return lib;
  }

  REACLIBReactionLibrary(const std::string& reaclibPath,
      const ReactionType type, const bool calclateInverses,
      const LeptonMode leptonMode, const std::string& description,
      const NuclideLibrary& nuclib, const NetworkOptions& opts,
      const bool doScreening, const bool freezeAtMinT9 = true);

  REACLIBReactionLibrary(const REACLIB& reaclib, const std::string& source,
      const ReactionType type, const bool calclateInverses,
      const LeptonMode leptonMode, const std::string& description,
      const NuclideLibrary& nuclib, const NetworkOptions& opts,
      const bool doScreening, const bool freezeAtMinT9 = true);

  REACLIBReactionLibrary(const ReactionType type, const bool calclateInverses,
      const std::string& description, const std::string& source,
      const std::vector<Reaction>& reactions,
      const std::vector<REACLIBCoefficients_t>& REACLIBCoefficients,
      const NuclideLibrary& nuclib, const NetworkOptions& opts,
      const bool doScreening, const bool freezeAtMinT9 = true);

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new REACLIBReactionLibrary(*this));
  }

  virtual void PrintAdditionalInfo(NetworkOutput * const pOutput) const;

  bool StateChanged(const ThermodynamicState thermoState,
      NetworkOutput * const pOutput);

  std::string Name() const {
    return "REACLIB Reaction Library";
  }

  std::size_t NumActiveNonConstReactions() const {
    return mNumActiveNonConstReactions;
  }

  bool CalcInverses() const {
    return mCalcInverses;
  }

  bool FreezeAtMinT9() const {
    return mFreezeAtMinT9;
  }

  const std::array<ReactionData<double>, NumREACLIBCoefficients>&
  RateFittingCoefficients() const {
    return mRateFittingCoefficients;
  }

  const ReactionData<bool>& RateIsTemperatureDependent() const {
    return mRateIsTemperatureDependent;
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

  template <bool SCR>
  void DoCalculateRatesScreening(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

private:
  REACLIBReactionLibrary(const NetworkOptions& opts):
      ReactionLibraryBase(opts) {}

  friend class boost::serialization::access;
  template <class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & boost::serialization::base_object<ReactionLibraryBase>(*this);
    ar & mNumActiveNonConstReactions;
    ar & mCalcInverses;
    ar & mFreezeAtMinT9;
    ar & mTempDepReactionsEnabled;
    for (unsigned int i = 0; i < NumREACLIBCoefficients; ++i)
      ar & mRateFittingCoefficients[i];
    ar & mRateIsTemperatureDependent;
    ar & mRates;
    ar & mInverseRates;
    ar & mPartfs;
    ar & mIsInverseReaction;
    ar & mGammas;
    ar & mDeltaNs;
    ar & mInverseRateSpinCorrections;
  }

  void CalculateREACLIBRateExpArguments(const double T9);

  void CalculateConstantRates();

  std::size_t mNumActiveNonConstReactions;

  // True if inverse rates should be calculated from detailed balance
  bool mCalcInverses;

  // True if temperature dependent rates should be frozen at the minimum
  // temperature for which they are valid, otherwise (if this is false), they
  // are turned off completely
  bool mFreezeAtMinT9;

  // state of the library (whether temperature dependent reactions are on or
  // off)
  bool mTempDepReactionsEnabled;

  // the REACLIB fitting coefficients of the rate
  std::array<ReactionData<double>, NumREACLIBCoefficients>
  mRateFittingCoefficients;

  ReactionData<bool> mRateIsTemperatureDependent;

  // the calculated reaction rates, the first mNumNonConstReactions values
  // are the non-constant calculated rates, the remaining values are the
  // constant rates
  std::vector<double> mRates;

  // inverse rates and partition function to correct inverse rates, these
  // are only used for strong reactions
  std::vector<double> mInverseRates;
  ReactionData<double> mPartfs;

  // only used if we are not calculating inverse rates ourselves
  ReactionData<bool> mIsInverseReaction;

  // TODO maybe compute these on the fly?
  // data for inverse reactions
  ReactionData<double> mGammas;
  ReactionData<int> mDeltaNs;
  ReactionData<double> mInverseRateSpinCorrections;
};


#endif // SKYNET_REACTIONS_REACLIBREACTIONLIBRARY_HPP_
