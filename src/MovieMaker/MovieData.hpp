/// \file MovieData.hpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_MOVIEDATA_HPP_
#define SKYNET_MOVIEMAKER_MOVIEDATA_HPP_

#ifdef SKYNET_USE_MOVIE

#include <limits>
#include <string>
#include <vector>

#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/Screening.hpp"
#include "Network/NetworkOutput.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

struct MinMax {
  MinMax():
    min(std::numeric_limits<int>::max()),
    max(0) { }

  int min, max;
};

class MovieData {
public:
  MovieData(const NetworkOutput& networkOutput, const double timeStep,
      const NeutrinoHistory* nuHist = nullptr,
      const Screening* screen = nullptr,
      const std::vector<const ReactionLibraryBase*>& reactionLibraries = {});

  MovieData(const NetworkOutput& networkOutput, const double initialTime,
      const double finalTime, const double timeStep,
      const NeutrinoHistory* nuHist = nullptr,
      const Screening* screen = nullptr,
      const std::vector<const ReactionLibraryBase*>& reactionLibraries = {});

  MovieData(const NetworkOutput& networkOutput, const double initialTime,
      const double finalTime, const int numTimes,
      const NeutrinoHistory* nuHist = nullptr,
      const Screening* screen = nullptr,
      const std::vector<const ReactionLibraryBase*>& reactionLibraries = {});

  const std::vector<int>& A() const {
    return mA;
  }

  const std::vector<std::vector<double>>&
  AbundanceVsAVsTLinInterp() const {
    return mAbundanceVsAVsTLinInterp;
  }

  const std::vector<std::vector<double>>&
  AbundanceVsIsotopeVsTLinInterp() const {
    return mAbundanceVsIsotopeVsTLinInterp;
  }

  const std::vector<std::vector<double>>&
  MassAbundanceVsAVsTLinInterp() const {
    return mMassAbundanceVsAVsTLinInterp;
  }

  const std::vector<std::vector<double>>&
  MassAbundanceVsIsotopeVsTLinInterp() const {
    return mMassAbundanceVsIsotopeVsTLinInterp;
  }

  const std::vector<std::vector<double>>&
  AbundanceVsAVsTLogInterp() const {
    return mAbundanceVsAVsTLogInterp;
  }

  const std::vector<std::vector<double>>&
  AbundanceVsIsotopeVsTLogInterp() const {
    return mAbundanceVsIsotopeVsTLogInterp;
  }

  const std::vector<std::vector<double>>&
  MassAbundanceVsAVsTLogInterp() const {
    return mMassAbundanceVsAVsTLogInterp;
  }

  const std::vector<std::vector<double>>&
  MassAbundanceVsIsotopeVsTLogInterp() const {
    return mMassAbundanceVsIsotopeVsTLogInterp;
  }

  const std::vector<std::vector<double>>&
  YDivTauBetaVsIsotopeVsTLogInterp() const {
    return mYDivTauBetaVsIsotopeVsTLogInterp;
  }

  const std::vector<double>& BindingEnergy() const {
    return mBindingEnergy;
  }

  const std::vector<double>& EdotVsT() const {
    return mEdotVsT;
  }

  const std::vector<double>& EntropyVsT() const {
    return mEntropyVsT;
  }

  int MaxA() const {
    return mMaxA;
  }

  int MaxZ() const {
    return mMaxZ;
  }

  int MaxN() const {
    return mMaxN;
  }

  int FirstTimeIndex() const {
    return mFirstTimeIndex;
  }

  int NumIsotopes() const {
    return mNumIsotopes;
  }

  int NumTimes() const {
    return mNumTimes;
  }

  const std::vector<double>& RhoVsT() const {
    return mRhoVsT;
  }

  const std::vector<double>& TemperatureVsT() const {
    return mTemperatureVsT;
  }

  const std::vector<double>& Time() const {
    return mTime;
  }

  const std::vector<double>& YeVsT() const {
    return mYeVsT;
  }

  const std::vector<int>& Z() const {
    return mZ;
  }

  const std::vector<int>& N() const {
    return mN;
  }

  const std::vector<MinMax>& ALimitsVsN() const {
    return mALimitsVsN;
  }

  const std::vector<MinMax>& ALimitsVsZ() const {
    return mALimitsVsZ;
  }

  const std::vector<MinMax>& NLimitsVsA() const {
    return mNLimitsVsA;
  }

  const std::vector<MinMax>& NLimitsVsZ() const {
    return mNLimitsVsZ;
  }

  const std::vector<MinMax>& ZLimitsVsA() const {
    return mZLimitsVsA;
  }

  const std::vector<MinMax>& ZLimitsVsN() const {
    return mZLimitsVsN;
  }

  bool HavePathInfo() const {
    return mHavePathInfo;
  }

private:
  void MakeVariousAbundances(
      const std::vector<std::vector<double>>& abundanceVsIsotopeVsT,
      std::vector<std::vector<double>> * const pMassAbundanceVsIsotopeVsT,
      std::vector<std::vector<double>> * const pAbundanceVsAVsT,
      std::vector<std::vector<double>> * const pMassAbundanceVsAVsT) const;

  int mNumIsotopes;
  int mNumTimes;
  int mFirstTimeIndex;
  int mMaxA;
  int mMaxZ;
  int mMaxN;
  bool mHavePathInfo;

  std::vector<int> mA;
  std::vector<int> mZ;
  std::vector<int> mN;
  std::vector<double> mBindingEnergy;
  std::vector<double> mTime;
  std::vector<double> mRhoVsT;
  std::vector<double> mTemperatureVsT;
  std::vector<double> mYeVsT;
  std::vector<double> mEntropyVsT;
  std::vector<double> mEdotVsT;
  std::vector<std::vector<double>> mAbundanceVsIsotopeVsTLinInterp;
  std::vector<std::vector<double>> mAbundanceVsAVsTLinInterp;
  std::vector<std::vector<double>> mMassAbundanceVsIsotopeVsTLinInterp;
  std::vector<std::vector<double>> mMassAbundanceVsAVsTLinInterp;
  std::vector<std::vector<double>> mAbundanceVsIsotopeVsTLogInterp;
  std::vector<std::vector<double>> mAbundanceVsAVsTLogInterp;
  std::vector<std::vector<double>> mMassAbundanceVsIsotopeVsTLogInterp;
  std::vector<std::vector<double>> mMassAbundanceVsAVsTLogInterp;
  std::vector<std::vector<double>> mYDivTauBetaVsIsotopeVsTLogInterp;

  std::vector<MinMax> mZLimitsVsN;
  std::vector<MinMax> mALimitsVsN;
  std::vector<MinMax> mNLimitsVsZ;
  std::vector<MinMax> mALimitsVsZ;
  std::vector<MinMax> mNLimitsVsA;
  std::vector<MinMax> mZLimitsVsA;
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_MOVIEDATA_HPP_
