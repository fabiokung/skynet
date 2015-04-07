/// \file ReactionPostProcess.cpp
/// \author lroberts
/// \since Oct 7, 2015
///
/// \brief
///
///

#include "Reactions/ReactionPostProcess.hpp"

ReactionPostProcess::ReactionPostProcess(
    const ReactionLibraryBase* reactionLib,
    const Screening * const pScreening,
    const NuclideLibrary& nucLib,
    const std::vector<ThermodynamicState>& thermoStates,
    const std::vector<std::vector<double>>& yIn,
    const std::vector<int>& ZIn,
    const std::vector<int>& AIn,
    const std::vector<double> times) :
    mpLibrary(reactionLib->MakeUniquePtr()),
    mpScreening(pScreening != nullptr ? pScreening->MakeUniquePtr() : nullptr),
    mTimes(times),
    mNuclib(nucLib) {
  if (mpScreening != nullptr)
    mpScreening->SetMaxZ(mpLibrary->MaxScreeningZ(ZIn));

  CalculateRates(thermoStates, yIn, ZIn, AIn);
}

ReactionPostProcess::ReactionPostProcess(
    const ReactionLibraryBase* reactionLib,
    const Screening * const pScreening,
    const NuclideLibrary& nucLib,
    const NetworkOutput& netOut,
    std::shared_ptr<NeutrinoHistory> nuHist) :
    mpLibrary(reactionLib->MakeUniquePtr()),
    mpScreening(pScreening != nullptr ? pScreening->MakeUniquePtr() : nullptr),
    mTimes(std::vector<double>()),
    mNuclib(nucLib) {
  std::vector<ThermodynamicState> thermoStates;

  for (unsigned int i = 0; i < netOut.Times().size(); ++i) {
    double time = netOut.Times()[i];
    mTimes.push_back(time);
    thermoStates.push_back(ThermodynamicState(
        netOut.TemperatureVsTime()[i],
        netOut.DensityVsTime()[i],
        netOut.EntropyVsTime()[i],
        std::numeric_limits<double>::signaling_NaN(), // U
        netOut.EtaEVsTime()[i],
        netOut.YeVsTime()[i],
        std::numeric_limits<double>::signaling_NaN(), // dSdT9
        (*nuHist)(time)));
  }

  if (mpScreening != nullptr)
    mpScreening->SetMaxZ(mpLibrary->MaxScreeningZ(netOut.Z()));

  CalculateRates(thermoStates, netOut.YVsTime(), netOut.Z(), netOut.A());
}

std::vector<std::vector<double>>
ReactionPostProcess::GetTotalNucleusWeakSingleReactantRate() const {
  std::vector<std::vector<const std::vector<double>*>> ratesVsTimeVsIsotope(
      mNuclib.NumNuclides());
  for (unsigned int r = 0; r < mpLibrary->NumActiveReactions(); ++r) {
    auto reaction = mpLibrary->Reactions()[r];

    if ((!reaction.IsWeak()) || (reaction.ReactantNames().size() != 1)
        || (reaction.NsOfReactants()[0] != 1)) {
      continue;
    }

    int nuclideId = mNuclib.NuclideIdsVsNames().at(reaction.ReactantNames()[0]);
    ratesVsTimeVsIsotope[nuclideId].push_back(&(mRates[r]));
  }

  std::vector<std::vector<double>> totalRatesVsIsotopeVsTime(
      mTimes.size());

  for (unsigned int i = 0; i < mTimes.size(); ++i)
    totalRatesVsIsotopeVsTime[i] = std::vector<double>(mNuclib.NumNuclides(),
        0.0);

  for (int j = 0; j < mNuclib.NumNuclides(); ++j) {
    for (auto rateVsTime : ratesVsTimeVsIsotope[j]) {
      for (unsigned int i = 0; i < mTimes.size(); ++i) {
        totalRatesVsIsotopeVsTime[i][j] += (*rateVsTime)[i];
      }
    }
  }

  return totalRatesVsIsotopeVsTime;
}

void ReactionPostProcess::CalculateRates(
    const std::vector<ThermodynamicState>& thermoStates,
    const std::vector<std::vector<double>>& yIn,
    const std::vector<int>& Z,
    const std::vector<int>& A) {

  assert(yIn[0].size() == (size_t)mNuclib.NumNuclides());
  assert(mTimes.size() == yIn.size());
  assert(mpLibrary->NumAllReactions() == mpLibrary->NumActiveReactions());

  std::map<std::string, int> specMap;
  for (unsigned int i = 0; i < Z.size(); i++) {
    specMap.insert(std::pair<std::string, int>(
        Nuclide::GetName(Z[i], A[i]), i));
  }

  mRates = std::vector<std::vector<double>>(mpLibrary->NumAllReactions(),
      std::vector<double>(mTimes.size()));

  mInvRates = std::vector<std::vector<double>>(mpLibrary->NumAllReactions(),
      std::vector<double>(mTimes.size()));

  // Calculate all reactions for each time
  for (unsigned int i = 0; i < mTimes.size(); ++i) {
    if (mpScreening != nullptr) {
      auto state = mpScreening->Compute(thermoStates[i], yIn[i], mNuclib);

      mpLibrary->CalculateRates(thermoStates[i],
          mNuclib.PartitionFunctionsWithoutSpinTerm(thermoStates[i].T9(), 100.0),
          100.0, &Z, &state.Mus());
    } else {
      mpLibrary->CalculateRates(thermoStates[i],
          mNuclib.PartitionFunctionsWithoutSpinTerm(thermoStates[i].T9(),
          100.0), 100.0, nullptr, nullptr);
    }

    for (unsigned int r = 0; r < mpLibrary->NumActiveReactions(); ++r) {
      auto reaction = mpLibrary->Reactions()[r];

      /// \todo Check that we have the double counting factors correct here
      double rate = mpLibrary->Rates()[r];
      for (auto& name : reaction.ReactantNames())
        rate *= yIn[i][specMap.at(name)];
      //rate *= yIn[i][mNuclib.NuclideIdsVsNames().at(name)];
      mRates[r][i] = rate;

      if (mpLibrary->InverseRates().size() > 0) {
        double invRate = mpLibrary->InverseRates()[r];
        for (auto& name : reaction.ProductNames())
          invRate *= yIn[i][specMap.at(name)];
        //invRate *= yIn[i][mNuclib.NuclideIdsVsNames().at(name)];
        mInvRates[r][i] = invRate;
      }
    }
  }
}

