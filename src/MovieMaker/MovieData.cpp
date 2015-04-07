/// \file MovieData.cpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#include "MovieMaker/MovieData.hpp"

#include <cassert>
#include <cmath>
#include <cstring>
#include <memory>
#include <valarray>

#include "Reactions/ReactionPostProcess.hpp"

namespace { // unnamed so this can only be used in this file

int GetNumTimes(const double initialTime, const double finalTime,
    const double deltaLogT) {
  double num = (log(finalTime) - log(initialTime)) / deltaLogT;
  return (int)(num + 1.5);
}

}

MovieData::MovieData(const NetworkOutput& networkOutput,
    const double initialTime, const double finalTime, const double timeStep,
    const NeutrinoHistory* nuHist, const Screening* screen,
    const std::vector<const ReactionLibraryBase*>& reactionLibraries) :
    MovieData(networkOutput, initialTime, finalTime,
        GetNumTimes(initialTime, finalTime, timeStep), nuHist,
        screen, reactionLibraries) {}

MovieData::MovieData(const NetworkOutput& networkOutput, const double timeStep,
    const NeutrinoHistory* nuHist, const Screening* screen,
    const std::vector<const ReactionLibraryBase*>& reactionLibraries):
    MovieData(networkOutput, networkOutput.Times()[0],
        networkOutput.Times()[networkOutput.NumEntries() - 1],
        GetNumTimes(networkOutput.Times()[0],
            networkOutput.Times()[networkOutput.NumEntries() - 1], timeStep),
            nuHist, screen, reactionLibraries) {}

MovieData::MovieData(const NetworkOutput& networkOutput,
    const double initialTime, const double finalTime, const int numTimes,
    const NeutrinoHistory* nuHist, const Screening* screen,
    const std::vector<const ReactionLibraryBase*>& reactionLibraries) {
  mA = networkOutput.A();
  mZ = networkOutput.Z();
  mBindingEnergy = networkOutput.BindingEnergy();

  mNumIsotopes = mA.size();
  assert((int)mZ.size() == mNumIsotopes);
  assert((int)mBindingEnergy.size() == mNumIsotopes);

  mN.resize(mNumIsotopes);
  for (int i = 0; i < mNumIsotopes; ++i)
    mN[i] = mA[i] - mZ[i];

  mMaxA = 0;
  mMaxZ = 0;
  mMaxN = 0;
  for (int i = 0; i < mNumIsotopes; ++i) {
    if (mA[i] > mMaxA)
      mMaxA = mA[i];
    if (mZ[i] > mMaxZ)
      mMaxZ = mZ[i];
    if (mN[i] > mMaxN)
      mMaxN = mN[i];
  }

  // get limits
  mZLimitsVsN.resize(mMaxN + 1);
  mALimitsVsN.resize(mMaxN + 1);
  mNLimitsVsZ.resize(mMaxZ + 1);
  mALimitsVsZ.resize(mMaxZ + 1);
  mNLimitsVsA.resize(mMaxA + 1);
  mZLimitsVsA.resize(mMaxA + 1);

  for (int i = 0; i < mNumIsotopes; ++i) {
    int n = (int)mN[i];
    int z = (int)mZ[i];
    int a = n + z;

    if (n > mNLimitsVsZ[z].max) mNLimitsVsZ[z].max = n;
    if (n < mNLimitsVsZ[z].min) mNLimitsVsZ[z].min = n;
    if (n > mNLimitsVsA[a].max) mNLimitsVsA[a].max = n;
    if (n < mNLimitsVsA[a].min) mNLimitsVsA[a].min = n;

    if (z > mZLimitsVsN[n].max) mZLimitsVsN[n].max = z;
    if (z < mZLimitsVsN[n].min) mZLimitsVsN[n].min = z;
    if (z > mZLimitsVsA[a].max) mZLimitsVsA[a].max = z;
    if (z < mZLimitsVsA[a].min) mZLimitsVsA[a].min = z;

    if (a > mALimitsVsN[n].max) mALimitsVsN[n].max = a;
    if (a < mALimitsVsN[n].min) mALimitsVsN[n].min = a;
    if (a > mALimitsVsZ[z].max) mALimitsVsZ[z].max = a;
    if (a < mALimitsVsZ[z].min) mALimitsVsZ[z].min = a;
  }

  // make times
  mFirstTimeIndex = 0;
  mTime = std::vector<double>(numTimes);
  double deltaLogT = (log(finalTime) - log(initialTime)) / (numTimes - 1);
  for (int i = 0; i < numTimes; ++i) {
    mTime[i] = exp(log(initialTime) + (double)i * deltaLogT);
  }

  auto inputTimes = networkOutput.Times();
  PiecewiseLinearFunction rhoFunc(inputTimes,
      networkOutput.DensityVsTime(), true);
  PiecewiseLinearFunction tempFunc(inputTimes,
      networkOutput.TemperatureVsTime(), true);
  PiecewiseLinearFunction yeFunc(inputTimes,
      networkOutput.YeVsTime(), true);
  PiecewiseLinearFunction entropyFunc(inputTimes,
      networkOutput.EntropyVsTime(), true);
  // don't do interpolation in log space because heating rate can be 0 or
  // negative
  PiecewiseLinearFunction edotFunc(inputTimes,
      networkOutput.HeatingRateVsTime(), false);

  mNumTimes = mTime.size();

  std::vector<std::valarray<double>> yVsTimeValarray(inputTimes.size());

  mHavePathInfo = (nuHist != nullptr) && (reactionLibraries.size() > 0);

  std::vector<std::vector<std::vector<double>>> tauBetas;
  std::vector<std::valarray<double>> yDivTauBetaVsTimeValarray(
      inputTimes.size());

  if (mHavePathInfo) {
    auto nuHistShrdPtr = nuHist->MakeSharedPtr();
    auto nuclib = networkOutput.GetNuclideLibrary();

    for (unsigned int i = 0; i < reactionLibraries.size(); ++i) {
      ReactionPostProcess reac(reactionLibraries[i], screen, nuclib,
          networkOutput, nuHistShrdPtr);
      tauBetas.push_back(reac.GetTotalNucleusWeakSingleReactantRate());
    }
  }

  for (unsigned int i = 0; i < inputTimes.size(); ++i) {
    std::valarray<double> valarray(networkOutput.YVsTime()[i].data(),
        mNumIsotopes);
    yVsTimeValarray[i] = valarray;
    yDivTauBetaVsTimeValarray[i] = valarray;

    if (mHavePathInfo) {
      std::valarray<double> tauBeta(0.0, mNumIsotopes);
      for (unsigned int j = 0; j < tauBetas.size(); ++j)
        tauBeta += std::valarray<double>(tauBetas[j][i].data(), mNumIsotopes);

      for (int j = 0; j < mNumIsotopes; ++j) {
        if ((tauBeta[j] >= 1.0E-50) && (valarray[j] >= 1.0E-18))
          yDivTauBetaVsTimeValarray[i][j] = valarray[j] * tauBeta[j];
        else
          yDivTauBetaVsTimeValarray[i][j] = 0.0;
      }
    }

    yVsTimeValarray[i] += 1.0E-100;
    yDivTauBetaVsTimeValarray[i] += 1.0E-100;
  }

  GeneralPiecewiseLinearFunction<std::valarray<double>> yFuncLin(inputTimes,
      yVsTimeValarray, false);
  GeneralPiecewiseLinearFunction<std::valarray<double>> yFuncLog(inputTimes,
      yVsTimeValarray, true);
  GeneralPiecewiseLinearFunction<std::valarray<double>> yDivTauBetaFuncLog(
      inputTimes, yDivTauBetaVsTimeValarray, true);

  mRhoVsT = std::vector<double>(mNumTimes);
  mTemperatureVsT = std::vector<double>(mNumTimes);
  mYeVsT = std::vector<double>(mNumTimes);
  mEntropyVsT = std::vector<double>(mNumTimes);
  mEdotVsT = std::vector<double>(mNumTimes);
  auto YVsTLin = std::vector<std::valarray<double>>(mNumTimes);
  auto YVsTLog = std::vector<std::valarray<double>>(mNumTimes);
  auto YDivTauBetaVsTLog = std::vector<std::valarray<double>>(mNumTimes);

  for (int i = 0; i < mNumTimes; ++i) {
    mRhoVsT[i] = rhoFunc(mTime[i]);
    mTemperatureVsT[i] = tempFunc(mTime[i]);
    mYeVsT[i] = yeFunc(mTime[i]);
    mEntropyVsT[i] = entropyFunc(mTime[i]);
    mEdotVsT[i] = edotFunc(mTime[i]);
    YVsTLin[i] = yFuncLin(mTime[i]);
    YVsTLog[i] = yFuncLog(mTime[i]);
    YDivTauBetaVsTLog[i] = yDivTauBetaFuncLog(mTime[i]);
  }

  mAbundanceVsIsotopeVsTLinInterp =
      std::vector<std::vector<double>>(mNumIsotopes);
  mAbundanceVsIsotopeVsTLogInterp =
      std::vector<std::vector<double>>(mNumIsotopes);
  mYDivTauBetaVsIsotopeVsTLogInterp =
      std::vector<std::vector<double>>(mNumIsotopes);

  for (int i = 0; i < mNumIsotopes; ++i) {
    mAbundanceVsIsotopeVsTLinInterp[i] = std::vector<double>(mNumTimes);
    mAbundanceVsIsotopeVsTLogInterp[i] = std::vector<double>(mNumTimes);
    mYDivTauBetaVsIsotopeVsTLogInterp[i] = std::vector<double>(mNumTimes);
  }

  for (int t = 0; t < mNumTimes; ++t) {
    for (int i = 0; i < mNumIsotopes; ++i) {
      mAbundanceVsIsotopeVsTLinInterp[i][t] = YVsTLin[t][i];
      mAbundanceVsIsotopeVsTLogInterp[i][t] = YVsTLog[t][i];
      mYDivTauBetaVsIsotopeVsTLogInterp[i][t] = YDivTauBetaVsTLog[t][i];
    }
  }

  MakeVariousAbundances(mAbundanceVsIsotopeVsTLinInterp,
      &mMassAbundanceVsIsotopeVsTLinInterp, &mAbundanceVsAVsTLinInterp,
      &mMassAbundanceVsAVsTLinInterp);

  MakeVariousAbundances(mAbundanceVsIsotopeVsTLogInterp,
      &mMassAbundanceVsIsotopeVsTLogInterp, &mAbundanceVsAVsTLogInterp,
      &mMassAbundanceVsAVsTLogInterp);
}

void MovieData::MakeVariousAbundances(
    const std::vector<std::vector<double>>& abundanceVsIsotopeVsT,
    std::vector<std::vector<double>> * const pMassAbundanceVsIsotopeVsT,
    std::vector<std::vector<double>> * const pAbundanceVsAVsT,
    std::vector<std::vector<double>> * const pMassAbundanceVsAVsT) const {
  // calculate mass abundances and abundaces vs A
  (*pMassAbundanceVsIsotopeVsT).resize(mNumIsotopes);
  (*pAbundanceVsAVsT).resize(mMaxA + 1);
  (*pMassAbundanceVsAVsT).resize(mMaxA + 1);

  for (int i = 0; i < mNumIsotopes; ++i) {
    assert((int)abundanceVsIsotopeVsT[i].size() == mNumTimes);
    (*pMassAbundanceVsIsotopeVsT)[i].resize(mNumTimes);

    for (int j = 0; j < mNumTimes; ++j) {
      (*pMassAbundanceVsIsotopeVsT)[i][j] =
          (double)mA[i] * abundanceVsIsotopeVsT[i][j];
    }
  }

  for (int i = 0; i < (mMaxA + 1); ++i) {
    (*pAbundanceVsAVsT)[i].resize(mNumTimes);
    (*pMassAbundanceVsAVsT)[i].resize(mNumTimes);
    memset((*pAbundanceVsAVsT)[i].data(), 0, mNumTimes * sizeof(double));
    memset((*pMassAbundanceVsAVsT)[i].data(), 0, mNumTimes * sizeof(double));
  }

  for (int i = 0; i < mNumIsotopes; ++i) {
    int A = int(mA[i]);
    for (int t = 0; t < mNumTimes; ++t) {
      (*pAbundanceVsAVsT)[A][t] += abundanceVsIsotopeVsT[i][t];
      (*pMassAbundanceVsAVsT)[A][t] += (*pMassAbundanceVsIsotopeVsT)[i][t];
    }
  }
}
