/// \file ReactionLibraryBase.hpp
/// \author jlippuner
/// \since Jul 18, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACTIONLIBRARYBASE_HPP_
#define SKYNET_REACTIONS_REACTIONLIBRARYBASE_HPP_

#include <functional>
#include <limits>
#include <set>
#include <vector>

#include <boost/serialization/access.hpp>

#include "EquationsOfState/EOS.hpp"
#include "MatrixSolver/Jacobian.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/NetworkOutput.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionData.hpp"
#include "Reactions/ReactionDataTypes.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name ReactionType in Python instead of the enum values
// being global constants

#ifdef SWIG
%rename(ReactionType) ReactionTypeStruct;
#endif // SWIG

struct ReactionTypeStruct {
  enum Value {
    Weak,
    Strong
  };
};

typedef ReactionTypeStruct::Value ReactionType;

class ReactionLibraryBase {
public:
  virtual ~ReactionLibraryBase() {}

  virtual std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const =0;

  virtual std::string Name() const =0;

  virtual const std::vector<double>& Rates() const =0;

  virtual const std::vector<double>& InverseRates() const =0;

  virtual const std::vector<double>& HeatingRates() const {
    return mEmptyDouble;
  }

  virtual const std::vector<double>& HeatingInverseRates() const {
    return mEmptyDouble;
  }

  void CalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

  void RemoveReactions(const std::vector<Reaction>& reactionsToBeRemoved);

  void RemoveReactions(const ReactionData<Reaction>& reactionsToBeRemoved);

  void SetState(const ThermodynamicState thermoState);

  virtual void SetTime(const double) {}

  virtual bool StateChanged(const ThermodynamicState thermoState,
      NetworkOutput * const pOutput) =0;

  void ReIndexNuclides(const NuclideLibrary& nuclib);

  double GetEDotExternal(const std::vector<double>& Y) const;

  void AddYdotContributions(const std::vector<double>& Y,
      std::vector<double> * const pYdot) const;

  void AddIndicesOfNonZeroJacobianEntries(
      std::set<std::array<int, 2>> * const indices) const;

  void AddJacobianContributions(const std::vector<double>& Y,
      Jacobian * const pJ) const;

  void PrintInfo(NetworkOutput * const pOutput) const;

  virtual void PrintAdditionalInfo(NetworkOutput * const /*pOutput*/) const {}

  void PrintRates(NetworkOutput * const pOutput) const;

  void SetDoOutput(bool out) {
    mDoScreenOutput = out;
  }

  bool GetDoOutput() const {
    return mDoScreenOutput;
  }

  std::size_t NumAllReactions() const {
    return mReactions.AllData().size();
  }

  std::size_t NumActiveReactions() const {
    return mReactions.ActiveData().size();
  }

  ReactionType Type() const {
    return mType;
  }

  const std::string& Description() const {
    return mDescription;
  }

  const std::string& Source() const {
    return mSource;
  }

  bool DoScreening() const {
    return mDoScreening;
  }

  const ReactionData<Reaction>& Reactions() const {
    return mReactions;
  }

  const ReactionData<ReactionData_t>& Data() const {
    return mData;
  }

  const ReactionData<ReactionNs4_t>& Ns0() const {
    return mNs0;
  }

  const ReactionData<ReactionNs4_t>& Ns1() const {
    return mNs1;
  }

  const ReactionData<ReactionIds4_t>& Ids0() const {
    return mIds0;
  }

  const ReactionData<ReactionIds4_t>& Ids1() const {
    return mIds1;
  }

  const std::vector<bool>& ActiveFlags() const {
    return mActiveFlags.AllData();
  }

  const ReactionData<double>& Qs() const {
    return mQs;
  }

  const ReactionData<std::string>& Labels() const {
    return mLabels;
  }

  int MaxScreeningZ(const std::vector<int>& Zs);

  const NetworkOptions& Options() const {
    return mOpts;
  }

protected:
  ReactionLibraryBase(const NetworkOptions& opts) :
      mOpts(opts) {}

  ReactionLibraryBase(const ReactionType type, const std::string& description,
      const std::string& source, const std::vector<Reaction>& reactions,
      const NuclideLibrary& nuclib, const NetworkOptions& opts,
      const bool doScreening);

  void UpdateActive(const std::vector<bool>& activeFlags);

  virtual void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func) =0;

  virtual void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const
          pScreeningChemicalPotentialCorrection) =0;

  const ReactionData<ReactionFactorials_t>& Factorials() const {
    return mFactorials;
  }

private:
  void LoopOverReactionDataInBase(
      const std::function<void(GeneralReactionData * const)>& func);

  void DoAddYdotContributions(const std::vector<double>& Y,
      std::vector<double> * const pYdot, const std::vector<double>& rates,
      const std::vector<double>& inverseRates) const;

  double DoGetEDotExternal(const std::vector<double>& Y,
      const std::vector<double>& rates,
      const std::vector<double>& inverseRates) const;

  void DoAddIndicesOfNonZeroJacobianEntries(
      std::set<std::array<int, 2>> * const indices,
      const bool alsoDoInverse) const;

  void DoAddJacobianContributions(const std::vector<double>& Y,
      Jacobian * const pJ, const std::vector<double>& rates,
      const std::vector<double>& inverseRates) const;

  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & mType;
    ar & mDescription;
    ar & mSource;
    ar & mStateHasBeenSet;
    ar & mDoScreening;
    ar & mReactions;
    ar & mActiveFlags;
    ar & mData;
    ar & mNs0;
    ar & mNs1;
    ar & mIds0;
    ar & mIds1;
    ar & mFactorials;
    ar & mQs;
    ar & mLabels;
    mDoScreenOutput = false;
  }

  // some meta data about the reactions stored in this library
  ReactionType mType;
  std::string mDescription;
  std::string mSource;

  bool mStateHasBeenSet;

  bool mDoScreening; // true if we apply screening corrections

  ReactionData<Reaction> mReactions;
  ReactionData<bool> mActiveFlags;

  // meta data of the reactions
  ReactionData<ReactionData_t> mData;

  // the numbers of reactants or products
  ReactionData<ReactionNs4_t> mNs0;
  ReactionData<ReactionNs4_t> mNs1;

  // the ids of the reactants or products
  ReactionData<ReactionIds4_t> mIds0;
  ReactionData<ReactionIds4_t> mIds1;

  // the products of factorials
  ReactionData<ReactionFactorials_t> mFactorials;

  // Q values
  ReactionData<double> mQs;

  // labels and descriptions
  ReactionData<std::string> mLabels;

  bool mDoScreenOutput;

  NetworkOptions mOpts;

  std::vector<double> mEmptyDouble;
};

#endif // SKYNET_REACTIONS_REACTIONLIBRARYBASE_HPP_
