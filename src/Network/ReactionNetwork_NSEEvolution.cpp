/// \file ReactionNetwork_NSEEvolution.cpp
/// \author jlippuner
/// \since Oct 31, 2014
///
/// \brief
///
///

#include "Network/ReactionNetwork.hpp"

#include "Network/NSE.hpp"

ReactionNetwork::YeEntropy& ReactionNetwork::YeEntropy::operator+=(
    const YeEntropy& other) {
  Ye += other.Ye;
  Entropy += other.Entropy;
  return *this;
}

ReactionNetwork::YeEntropy ReactionNetwork::YeEntropy::operator*(
    const double a) const {
  YeEntropy res;
  res.Ye = Ye * a;
  res.Entropy = Entropy * a;
  return res;
}

void ReactionNetwork::SetupNSEEvolution() {
  // set up stuff for NSE evolution
  if (mOpts.IsSelfHeating) {
    mpNseEvolutionRKF45 = std::unique_ptr<RKF45<YeEntropy>>(
      new RKF45<YeEntropy>([&] (const double time, const YeEntropy& val) {
        double rho = (*mpRhoVsTime)(time);
        auto state = mCurrentThermodynamicState;

        mProfiler.NSEEvolution.Start();
        auto nseResult = mpNSE->CalcFromEntropyAndDensityWithTemperatureGuess(
            val.Entropy, rho, val.Ye, state.T9());
        mProfiler.NSEEvolution.Stop();
        state = mpEOS->FromTAndRho(nseResult.T9(), rho, nseResult.Y(),
            *mpNuclideLibrary, (*mpNuDistVsTime)(time));

        this->CalculateRates(time, state);
        std::vector<double> yDot(nseResult.Y().size(), 0.0);

        this->AddWeakYdotContributions(time, nseResult.Y(), &yDot);

        double YeDot = 0.0;
        for (unsigned int i = 0; i < yDot.size(); ++i)
          YeDot += mpNuclideLibrary->Zs()[i] * yDot[i];

        auto heatEntropy = CalculateHeatingRateAndEntropyChangeOrDerivative(
            time, 0.0, state.T9(), rho, mCurrentY, yDot, true);

        YeEntropy res;
        res.Ye = YeDot;
        res.Entropy = heatEntropy.EntropyChangeOrDerivative;

        return res;
      }));
  } else {
    mpNseEvolutionRKF45 = std::unique_ptr<RKF45<YeEntropy>>(
      new RKF45<YeEntropy>([&] (const double time, const YeEntropy& val) {
        double rho = (*mpRhoVsTime)(time);
        double t9 = (*mpT9VsTime)(time);
        auto state = mCurrentThermodynamicState;

        mProfiler.NSEEvolution.Start();
        auto nseResult = mpNSE->CalcFromTemperatureAndDensity(t9, rho, val.Ye);
        mProfiler.NSEEvolution.Stop();
        state = mpEOS->FromTAndRho(t9, rho, nseResult.Y(),
            *mpNuclideLibrary, (*mpNuDistVsTime)(time));

        this->CalculateRates(time, state);
        std::vector<double> yDot(nseResult.Y().size(), 0.0);

        this->AddWeakYdotContributions(time, nseResult.Y(), &yDot);

        double YeDot = 0.0;
        for (unsigned int i = 0; i < yDot.size(); ++i)
          YeDot += mpNuclideLibrary->Zs()[i] * yDot[i];

        YeEntropy res;
        res.Ye = YeDot;
        res.Entropy = 0.0;

        return res;
    }));
  }

  // TODO add options for NSE evolution step size
  if (mOpts.IsSelfHeating) {
    mpNseEvolutionIntegrator =
        std::unique_ptr<AdaptiveODEIntegrator<std::array<YeEntropy, 2>>>(
        new AdaptiveODEIntegrator<std::array<YeEntropy, 2>>(
            [&] (const std::array<YeEntropy, 2>& prev, const double time,
                const double dt) {
              return mpNseEvolutionRKF45->NextStep(prev[0], time, dt);
            },
            [&] (const std::array<YeEntropy, 2>& dat) {
              double Ye = dat[0].Ye;
              if ((Ye <= 0.0) || (Ye >= 1.0))
                return std::numeric_limits<double>::max();

              double YeError = 2.0 * fabs(dat[0].Ye - dat[1].Ye)
                  / fabs(dat[0].Ye + dat[1].Ye);

              double s = dat[0].Entropy;
              if (s <= 0.0)
                return std::numeric_limits<double>::max();

              double sError = 2.0 * fabs(dat[0].Entropy - dat[1].Entropy)
                  / fabs(dat[0].Entropy + dat[1].Entropy);

              return std::max(YeError, sError);
            },
            mOpts.MinDt, mOpts.MaxDt, mOpts.NSEEvolutionMinError,
            mOpts.NSEEvolutionMaxError, mOpts.MaxDtChangeMultiplier)
        );
  } else {
    mpNseEvolutionIntegrator =
    std::unique_ptr<AdaptiveODEIntegrator<std::array<YeEntropy, 2>>>(
    new AdaptiveODEIntegrator<std::array<YeEntropy, 2>>(
        [&] (const std::array<YeEntropy, 2>& prev, const double time,
            const double dt) {
          return mpNseEvolutionRKF45->NextStep(prev[0], time, dt);
        },
        [&] (const std::array<YeEntropy, 2>& dat) {
          double Ye = dat[0].Ye;
          if ((Ye <= 0.0) || (Ye >= 1.0))
            return std::numeric_limits<double>::max();

          double YeError = 2.0 * fabs(dat[0].Ye - dat[1].Ye)
              / fabs(dat[0].Ye + dat[1].Ye);

          return YeError;
        },
        mOpts.MinDt, mOpts.MaxDt, mOpts.NSEEvolutionMinError,
        mOpts.NSEEvolutionMaxError, mOpts.MaxDtChangeMultiplier)
    );
  }
}

double ReactionNetwork::CheckIfNSEEvolution() {
  // determine if we turn on or off NSE evolution
  if (mDoNSEEvolution
      && (mCurrentThermodynamicState.T9() < mOpts.NSEEvolutionMinT9)) {
    mDoNSEEvolution = false;
    mpOutput->Log("# Switching to full network evolution (temperature below "
        "threshold for NSE evolution)\n");
  }

  // find density time scale
  double densityTimeScaleInverse = std::numeric_limits<double>::max();
  if (mPreviousDtAndLimitString.Dt > 0.0) {
    double rho0 = (*mpRhoVsTime)(mCurrentTime);
    double rho1 = (*mpRhoVsTime)(mCurrentTime + mPreviousDtAndLimitString.Dt);
    densityTimeScaleInverse =
        fabs((rho1 - rho0) / (rho0 * mPreviousDtAndLimitString.Dt));
  }

  if (!mDoNSEEvolution
      && (mCurrentThermodynamicState.T9() >= mOpts.NSEEvolutionMinT9)) {
    // find longest strong interaction time scale
    double strongTimeScaleInverse = std::numeric_limits<double>::max();

    auto thermoStateForRates = mpEOS->FromTAndRho(
        mCurrentThermodynamicState.T9(), (*mpRhoVsTime)(mCurrentTime),
        mCurrentY, *mpNuclideLibrary, (*mpNuDistVsTime)(mCurrentTime));

    std::vector<double> yDotStrong(mCurrentY.size(), 0.0);
    auto partf = mpNuclideLibrary->PartitionFunctionsWithoutSpinTerm(
        thermoStateForRates.T9(), mOpts.PartitionFunctionLog10Cap);

    if (mpScreening != nullptr) {
      auto screen = mpScreening->Compute(thermoStateForRates, mCurrentY,
          *mpNuclideLibrary);

      for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i) {
        if (mpReactionLibraries[i]->Type() == ReactionType::Strong)
          mpReactionLibraries[i]->CalculateRates(thermoStateForRates, partf,
              mOpts.RateExpArgumentCap, &mpNuclideLibrary->Zs(), &screen.Mus());
        mpReactionLibraries[i]->AddYdotContributions(mCurrentY, &yDotStrong);
      }
    } else {
      for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i) {
        if (mpReactionLibraries[i]->Type() == ReactionType::Strong)
          mpReactionLibraries[i]->CalculateRates(thermoStateForRates, partf,
              mOpts.RateExpArgumentCap, nullptr, nullptr);
        mpReactionLibraries[i]->AddYdotContributions(mCurrentY, &yDotStrong);
      }
    }

    double smallestYdotByY = std::numeric_limits<double>::max();
#if defined(__ICC) || defined(__INTEL_COMPILER)
    // if the Intel compiler vectorizes this loop, it can trigger a
    // floating point exception when it divides by zero (even though
    // it would not use the result due to the if statement, the
    // floating point exception still triggers)
#pragma novector
#endif // Intel
    for (unsigned int i = 0; i < mCurrentY.size(); ++i) {
      if (mCurrentY[i] >= mOpts.SmallestYUsedForDtCalculation) {
        smallestYdotByY = std::min(smallestYdotByY,
            fabs(yDotStrong[i] / mCurrentY[i]));
      }
    }

    strongTimeScaleInverse = smallestYdotByY;

    if (strongTimeScaleInverse > mOpts.NSEEvolutionStartTimeScaleFactor
        * densityTimeScaleInverse) {
      SwitchToNSEEvolutionIfPossible(
          mOpts.MaxFractionalTemperatureAndEntropyChange, "strong interaction "
              "time scale is smaller than density time scale", nullptr);
    }
  }

  if (!mDoNSEEvolution &&
      mCurrentThermodynamicState.T9() >= mOpts.NetworkEvolutionMaxT9) {
    SwitchToNSEEvolutionIfPossible(
        mOpts.MaxFractionalTemperatureAndEntropyChange, "strong interaction "
            "time scale is smaller than density time scale", nullptr);
  }

  double dtMax = mOpts.MaxDt;
  if (densityTimeScaleInverse > 0.0)
    dtMax = mOpts.DensityTimeScaleMultiplier / densityTimeScaleInverse;

  return dtMax;
}

bool ReactionNetwork::SwitchToNSEEvolutionIfPossible(
    const double fractionalJumpLimit, const std::string& reason,
    std::vector<double> * const pY) {
  if (mCurrentThermodynamicState.T9() < mOpts.NSEEvolutionMinT9)
    return false;

  // calc NSE to see if we're close
  auto nse = mpNSE->CalcFromInternalEnergyAndDensityWithTemperatureGuess(
      mCurrentThermodynamicState.U(), mCurrentThermodynamicState.Rho(),
      mCurrentThermodynamicState.Ye(), mCurrentThermodynamicState.T9());
  double nseS = mpEOS->FromTAndRho(nse.T9(), nse.Rho(), nse.Y(),
      *mpNuclideLibrary).S();

  // switch to NSE evolution if temperatures are close
  if ((fabs(nse.T9() - mCurrentThermodynamicState.T9())
      / mCurrentThermodynamicState.T9() > fractionalJumpLimit)
      || (fabs(nseS - mCurrentThermodynamicState.S())
          / mCurrentThermodynamicState.S() > fractionalJumpLimit)) {
    return false;
  } else {
    mDoNSEEvolution = true;
    mpOutput->Log("# Switching to NSE evolution (%s)\n", reason.c_str());

    // set initial nseStep
    mNseStep.t = mCurrentTime;
    mNseStep.nextDt = mPreviousDtAndLimitString.Dt;
    mNseStep.value[0].Ye = mCurrentThermodynamicState.Ye();
    double nseEntropy = mpEOS->FromTAndRho(nse.T9(), nse.Rho(), nse.Y(),
        *mpNuclideLibrary).S();
    mNseStep.value[0].Entropy = nseEntropy;
    mNseStep.value[1] = mNseStep.value[0];

    if (pY != nullptr)
      *pY = nse.Y();

    return true;
  }
}

void ReactionNetwork::EvolveNSE(const double dtMax) {
  mNseStep = mpNseEvolutionIntegrator->TakeStep(mNseStep, dtMax);
  double Ye = mNseStep.value[0].Ye;

  double rho = (*mpRhoVsTime)(mNseStep.t);
  mProfiler.NSEEvolution.Start();
  std::unique_ptr<NSEResult> pNseRes;
  if (mOpts.IsSelfHeating) {
    pNseRes = std::unique_ptr<NSEResult>(
        new NSEResult(mpNSE->CalcFromEntropyAndDensityWithTemperatureGuess(
            mNseStep.value[0].Entropy, rho, Ye,
            mCurrentThermodynamicState.T9())));
  } else {
    double t9 = (*mpT9VsTime)(mNseStep.t);
    pNseRes = std::unique_ptr<NSEResult>(
        new NSEResult(mpNSE->CalcFromTemperatureAndDensity(t9, rho, Ye)));
  }
  mProfiler.NSEEvolution.Stop();

  mPreviousY.swap(mCurrentY);
  mCurrentY = pNseRes->Y();
  mCurrentThermodynamicState = mpEOS->FromTAndRho(pNseRes->T9(), rho,
      mCurrentY, *mpNuclideLibrary, (*mpNuDistVsTime)(mNseStep.t) );

  std::vector<double> yDot(mCurrentY.size(), 0.0);
  AddWeakYdotContributions(mNseStep.t, mCurrentY, &yDot);

  auto heatEntropy = CalculateHeatingRateAndEntropyChangeOrDerivative(
      mNseStep.t, 0.0, mCurrentThermodynamicState.T9(),
      mCurrentThermodynamicState.Rho(), mCurrentY, yDot, true);
  mCurrentHeatingRate = heatEntropy.HeatingRate;

  mCurrentNumberOfFailedTimeSteps = 0;
  mCurrentNumberOfNRIterations = 0;
  mPreviousDtAndLimitString.Dt = mNseStep.t - mCurrentTime;
  mPreviousDtAndLimitString.Str = "nse";
}

