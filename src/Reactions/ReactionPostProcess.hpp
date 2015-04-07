/// \file ReactionPostProcess.hpp
/// \author lroberts
/// \since Oct 7, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACTIONPOSTPROCESS_HPP_
#define SKYNET_REACTIONS_REACTIONPOSTPROCESS_HPP_

#include <vector>
#include <string>

#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/Screening.hpp"
#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "NuclearData/Nuclide.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

struct RateVsTime {
  Reaction reac;
  std::vector<double> rate;
};

typedef const std::vector<std::string>&(Reaction::*NameFunc)() const;
typedef const std::vector<int>&(Reaction::*NumberFunc)() const;

class ReactionPostProcess {
public:
  ReactionPostProcess(
      const ReactionLibraryBase* reactionLib,
      const Screening * const pScreening,
      const NuclideLibrary& nucLib,
      const std::vector<ThermodynamicState>& thermoStates,
      const std::vector<std::vector<double>>& yIn,
      const std::vector<int>& ZIn,
      const std::vector<int>& AIn,
      const std::vector<double> times);

  ReactionPostProcess(
      const ReactionLibraryBase* reactionLib,
      const Screening * const pScreening,
      const NuclideLibrary& nucLib,
      const NetworkOutput& netOut,
      std::shared_ptr<NeutrinoHistory> nuHist);

  std::vector<double> GetTimes() {
    return mTimes;
  }

  std::vector<double> GetTotalNucleusDestructionRate(std::string nuc) {
    return GetTotalRate<&Reaction::ReactantNames, &Reaction::NsOfReactants,
        &Reaction::ProductNames, &Reaction::NsOfProducts>(nuc);
  }

  std::vector<double> GetTotalNucleusProductionRate(std::string nuc) {
    return GetTotalRate<&Reaction::ProductNames, &Reaction::NsOfProducts,
        &Reaction::ReactantNames, &Reaction::NsOfReactants>(nuc);
  }

  std::vector<std::vector<double>>
  GetTotalNucleusWeakSingleReactantRate() const;

  std::vector<RateVsTime> GetNucleusDestructionRates(std::string nuc) {
    return GetRates<&Reaction::ReactantNames, &Reaction::NsOfReactants,
        &Reaction::ProductNames, &Reaction::NsOfProducts>(nuc);
  }

  std::vector<RateVsTime> GetNucleusProductionRates(std::string nuc) {
    return GetRates<&Reaction::ProductNames, &Reaction::NsOfProducts,
        &Reaction::ReactantNames, &Reaction::NsOfReactants>(nuc);
  }

protected:
  template<NameFunc Names, NumberFunc Numbers, NameFunc NamesInv,
      NumberFunc NumberInv>
  std::vector<double> GetTotalRate(std::string nuc);

  template<NameFunc Names, NumberFunc Numbers, NameFunc NamesInv,
      NumberFunc NumberInv>
  std::vector<RateVsTime> GetRates(std::string nuc);

  void CalculateRates(const std::vector<ThermodynamicState>& thermoStates,
      const std::vector<std::vector<double>>& yIn,
      const std::vector<int>& Z,
      const std::vector<int>& A);

  std::unique_ptr<ReactionLibraryBase> mpLibrary;
  std::unique_ptr<Screening> mpScreening;
  std::vector<double> mTimes;
  NuclideLibrary mNuclib;
  std::vector<std::vector<double>> mRates, mInvRates;
};

// Have to define template functions in the header
template<NameFunc Names, NumberFunc Numbers,
    NameFunc NamesInv, NumberFunc NumbersInv>
std::vector<double> ReactionPostProcess::GetTotalRate(std::string nuc) {
  auto rates = GetRates<Names, Numbers, NamesInv, NumbersInv>(nuc);
  std::vector<double> totRate(mTimes.size(), 0.0);
  for (auto& rate : rates) {
    for (unsigned int i = 0; i < mTimes.size(); ++i) {
      totRate[i] += rate.rate[i];
    }
  }
  return totRate;
}

template<NameFunc Names, NumberFunc Numbers,
    NameFunc NamesInv, NumberFunc NumbersInv>
std::vector<RateVsTime> ReactionPostProcess::GetRates(std::string nuc) {
  std::vector<RateVsTime> rates;
  for (unsigned int r = 0; r < mpLibrary->NumActiveReactions(); ++r) {
    auto reaction = mpLibrary->Reactions()[r];

    double nreac = 0.0;
    for (unsigned int i = 0; i < (reaction.*Names)().size(); ++i) {
      if ((reaction.*Names)()[i] == nuc) {
        nreac = (double)(reaction.*Numbers)()[i];
        break;
      }
    }

    if (nreac > 0.0) {
      RateVsTime rvt;
      rvt.reac = reaction;
      rvt.rate = mRates[r];
      rates.push_back(rvt);
    }

    double nreacInv = 0.0;
    for (unsigned int i = 0; i < (reaction.*NamesInv)().size(); ++i) {
      if ((reaction.*NamesInv)()[i] == nuc) {
        nreacInv = (double)(reaction.*NumbersInv)()[i];
        break;
      }
    }

    if (nreacInv > 0.0) {
      RateVsTime rvt;
      rvt.reac = reaction;
      rvt.rate = mInvRates[r];
      rates.push_back(rvt);
    }

  }

  return rates;
}

#endif // SKYNET_REACTIONS_REACTIONPOSTPROCESS_HPP_
