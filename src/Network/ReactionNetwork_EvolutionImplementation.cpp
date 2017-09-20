/// \file ReactionNetwork.cpp
/// \author jlippuner
/// \since Jun 30, 2014
///
/// \brief
///
///

#include "Network/ReactionNetwork.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "DensityProfiles/ConstantFunction.hpp"
#include "MatrixSolver/MatrixSolver.hpp"
#include "Network/NSE.hpp"
#include "Utilities/Constants.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"

void ReactionNetwork::InitOperatorSplitEvolution(
    const std::vector<double>& initialY, const double startTime,
    const double initialT9, const double initialRho,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  ConstantFunction T9VsTime(initialT9);
  ConstantFunction rhoVsTime(initialRho);

  SetupEvolution(InitMode::TemperatureVsTime, initialY, startTime, NaN, initDt,
      NaN, NaN, &T9VsTime, &rhoVsTime, NaN, outputPrefix,
      massFractionOutputNuclideNames);
}

OperatorSplitStepResult ReactionNetwork::OperatorSplitTakeStep(const double T9,
    const double rho, const double endTime) {
  // set constant temperature and rho
  ConstantFunction T9VsTime(T9);
  ConstantFunction rhoVsTime(rho);
  mpT9VsTime = T9VsTime.MakeUniquePtr();
  mpRhoVsTime = rhoVsTime.MakeUniquePtr();

  // SET THE NEUTRINO DISTRIBUTION FUNCTION IF DESIRED //

  // save current time step index
  std::size_t lastTimeStepIdx = mpOutput->NumEntries();

  // take steps until endTime
  DoEvolutionSteps(endTime);

  // add up heating contribution and entropy changes from each step taken
  double totalHeating = 0.0; // in erg / gram
  for (unsigned int i = lastTimeStepIdx; i < mpOutput->NumEntries(); ++i)
    totalHeating += mpOutput->HeatingRateVsTime()[i] * mpOutput->DtVsTime()[i];

  // in k_B / baryon
  double totalEntropyChange =
      mpOutput->EntropyVsTime()[mpOutput->NumEntries() - 1]
          - mpOutput->EntropyVsTime()[lastTimeStepIdx - 1];

  return OperatorSplitStepResult(totalHeating, totalEntropyChange);
}

void ReactionNetwork::FinishOperatorSplitEvolution() {
  FinishEvolution();
}

const NetworkOutput& ReactionNetwork::DoEvolution(const InitMode initMode,
    const std::vector<double>& initialY,
    const double startTime, const double endTime, const double initDt,
    const double initialT9, const double initialEntropy,
    const FunctionVsTime<double> * const pT9VsTime,
    const FunctionVsTime<double> * const pRhoVsTime,
    const double Ye, const std::string& outputPrefix,
    std::vector<std::string> massFractionOutputNuclideNames) {
  SetupEvolution(initMode, initialY, startTime, endTime, initDt, initialT9,
      initialEntropy, pT9VsTime, pRhoVsTime, Ye, outputPrefix,
      massFractionOutputNuclideNames);
  DoEvolutionSteps(endTime);
  FinishEvolution();
  return *mpOutput;
}

void ReactionNetwork::SetupEvolution(const InitMode initMode,
    const std::vector<double>& initialY,
    const double startTime, const double endTime, const double initDt,
    const double initialT9, const double initialEntropy,
    const FunctionVsTime<double> * const pT9VsTime,
    const FunctionVsTime<double> * const pRhoVsTime,
    const double Ye, const std::string& outputPrefix,
    std::vector<std::string> massFractionOutputNuclideNames) {
  // check that network if self-heating if we don't have temperature vs. time
  if (initMode != InitMode::TemperatureVsTime) {
    if (!mOpts.IsSelfHeating) {
      throw std::runtime_error("Must specify a temperature vs. time history "
          "for a non-self-heating network");
    }
  }

  if (!((!std::isnan(Ye)) ^ (initialY.size() != 0))) {
    throw std::runtime_error("Must either set Ye to a number or provide "
        "initial abundances");
  }

  if (pRhoVsTime == nullptr) {
    throw std::invalid_argument("Pointer to rhoVsTime function must not be "
        "null");
  }

  // reset output
  if (massFractionOutputNuclideNames.size() == 0) {
    mpOutput = make_unique<NetworkOutput>(
        NetworkOutput::CreateNew(outputPrefix, *mpNuclideLibrary,
            mOpts.DisableStdoutOutput));
  } else {
    std::vector<int> massFractionOutputNuclideIds;
    for (auto& name : massFractionOutputNuclideNames) {
      massFractionOutputNuclideIds.push_back(
          mpNuclideLibrary->NuclideIdsVsNames().at(name));
    }

    mpOutput = make_unique<NetworkOutput>(
        NetworkOutput::CreateNew(outputPrefix, *mpNuclideLibrary,
            mOpts.DisableStdoutOutput, massFractionOutputNuclideIds));
  }

  int numWeakReactions = 0;
  int numStrongReactions = 0;
  int numInverseReactions = 0;

  for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i) {
    switch (mpReactionLibraries[i]->Type()) {
    case ReactionType::Weak:
      numWeakReactions += mpReactionLibraries[i]->NumActiveReactions();
      break;
    case ReactionType::Strong:
      numStrongReactions += mpReactionLibraries[i]->NumActiveReactions();
      break;
    default:
      throw std::runtime_error("Unknown reaction type");
    }

    numInverseReactions += mpReactionLibraries[i]->InverseRates().size();
  }

  mpOutput->Log("# Network contains %i nuclides and %i reactions "
      "(%i weak, %i strong, %i inverse)\n#\n", mpNuclideLibrary->NumNuclides(),
      numWeakReactions + numStrongReactions + numInverseReactions,
      numWeakReactions, numStrongReactions, numInverseReactions);

  mpNuclideLibrary->PrintInfo(mpOutput.get());

  for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i)
    mpReactionLibraries[i]->PrintInfo(mpOutput.get());

  mpEOS->PrintInfo(mpOutput.get());
  mpOutput->Log("#\n");

  if (mpScreening.get() == nullptr)
    mpOutput->Log("# Screening: no screening\n");
  else
    mpScreening->PrintInfo(mpOutput.get());
  mpOutput->Log("#\n");

  mpNuDistVsTime->PrintInfo(mpOutput.get());
  mpOutput->Log("#\n");

  mpJacobian->Describe(mpOutput.get());
  mpOutput->Log("#\n");

  mpOutput->Log("# Options:\n");
  mOpts.WriteToLog(mpOutput.get());
  mpOutput->Log("#\n");

  mpOutput->WriteColumnDescriptionToLog();

  // set initial time and density
  mCurrentTime = startTime;
  mpRhoVsTime = pRhoVsTime->MakeUniquePtr();
  double rho = (*mpRhoVsTime)(mCurrentTime);

  auto initialYToBeUsed = initialY;
  bool startFromNSE = false;

  if (!std::isnan(Ye)) {
    startFromNSE = true;
    // we need to do NSE
    if ((initMode == InitMode::TemperatureVsTime)
        || (initMode == InitMode::Temperature)) {
      // we have an initial temperature, so we can run NSE directly
      double T9 =
          (initMode == InitMode::TemperatureVsTime) ?
              (*pT9VsTime)(mCurrentTime) : initialT9;

      mpOutput->Log("# Calculating NSE with T = %.3E GK, rho = %.3E g / cm^3, "
          "Ye = %.3f\n", T9, rho, Ye);
      auto nseOut = mpNSE->CalcFromTemperatureAndDensity(T9, rho, Ye);
      mpOutput->Log("# NSE converged in %i iterations with an error of %.3E\n",
          nseOut.NumIterations(), nseOut.Error());
      mpOutput->Log("#\n");

      initialYToBeUsed = nseOut.Y();
    } else {
      mpOutput->Log("# Calculating NSE with entropy = %.3E k_B / baryon, "
          "rho = %.3E g / cm^3, Ye = %.3f\n", initialEntropy, rho, Ye);
      auto nseOut = mpNSE->CalcFromEntropyAndDensity(initialEntropy, rho, Ye);
      mpOutput->Log("# NSE found T = %.3E GK and converged in %i iterations "
          "with an error of %.3E\n", nseOut.T9(), nseOut.NumIterations(),
          nseOut.Error());
      mpOutput->Log("#\n");

      initialYToBeUsed = nseOut.Y();
    }
  }

  // set initial Y
  mCurrentY = std::vector<double>(initialYToBeUsed);
  mPreviousY = std::vector<double>(0);

  // renormalize initial Y so that sum of mass fractions is 1
  RenormalizeMass();

  // set initial temperature and entropy
  if ((initMode == InitMode::TemperatureVsTime)
      || (initMode == InitMode::Temperature)) {
    double T9ToUse = 0.0;
    if (initMode == InitMode::TemperatureVsTime) {
      mpT9VsTime = pT9VsTime->MakeUniquePtr();
      T9ToUse = (*mpT9VsTime)(mCurrentTime);
    } else {
      T9ToUse = initialT9;
    }

    mCurrentThermodynamicState = mpEOS->FromTAndRho(T9ToUse, rho, mCurrentY,
        *mpNuclideLibrary, (*mpNuDistVsTime)(mCurrentTime));
  } else {
    if (mOpts.IsSelfHeating)
      mCurrentThermodynamicState = mpEOS->FromSAndRho(initialEntropy, rho,
          mCurrentY, *mpNuclideLibrary, (*mpNuDistVsTime)(mCurrentTime));
    else
      throw std::runtime_error("Cannot set initial temperature from entropy "
          "for a non-self-heating network");
  }

  mPreviousDtAndLimitString = DtStr(0.0, "none");

  if (!std::isnan(endTime)) {
    mpOutput->LogWithTimestamp("# Starting evolution, end time = %.3E\n#\n",
        endTime);
  } else {
    // this is an operator split evolution
    mpOutput->LogWithTimestamp("# Starting operator split evolution\n#\n");
  }

  mCurrentHeatingRate = 0.0;

  // set initial state on the reaction libraries
  for (auto& lib : mpReactionLibraries)
    lib->SetState(mCurrentThermodynamicState);
  PostStepUpdate();

  // get initial dt
  if (!std::isnan(endTime)) {
    mPreviousDtAndLimitString =
        CalculateFirstDt(endTime - mCurrentTime, initDt);
  } else {
    // this is an operator split evolution
    mPreviousDtAndLimitString =
        CalculateFirstDt(mOpts.MaxDt, initDt);
  }

  // set up NSE evolution
  SetupNSEEvolution();

  if (mCurrentThermodynamicState.T9() >= mOpts.NetworkEvolutionMaxT9) {
    startFromNSE = true;
  }

  mDoNSEEvolution = (startFromNSE &&
      (mCurrentThermodynamicState.T9() >= mOpts.NSEEvolutionMinT9));
  if (mDoNSEEvolution) {
    // set initial nseStep
    mNseStep.t = mCurrentTime;
    mNseStep.nextDt = mPreviousDtAndLimitString.Dt;
    mNseStep.value[0].Ye = mCurrentThermodynamicState.Ye();
    mNseStep.value[0].Entropy = mCurrentThermodynamicState.S();
    mNseStep.value[1] = mNseStep.value[0];
  }

  mProfiler.Everything.Start();
}

void ReactionNetwork::DoEvolutionSteps(const double endTime) {
  while (mCurrentTime < endTime) {
    double nseEvolutionDtMax = CheckIfNSEEvolution();

    if (mDoNSEEvolution) {
      // NSE evolution step
      EvolveNSE(std::min(nseEvolutionDtMax, endTime - mCurrentTime));
    } else {
      // ordinary network evolution step
      mPreviousDtAndLimitString = CalculateDt(endTime - mCurrentTime);
      mPreviousDtAndLimitString.Dt = TakeStep(mPreviousDtAndLimitString.Dt,
          false);
    }

    mCurrentTime += mPreviousDtAndLimitString.Dt;
    PostStepUpdate();
  }
}

void ReactionNetwork::FinishEvolution() {
  mProfiler.Everything.Stop();
  mpOutput->Log("#\n");
  mpOutput->LogWithTimestamp("# Finished evolution\n#\n");

  mProfiler.Report(mpOutput.get(), mTotalNumberOfNRIterations);

  mpOutput->Close();
}

void ReactionNetwork::CalculateRates(const double time,
    const ThermodynamicState thermoState) {
  mProfiler.RateCalculation.Start();

  auto partf = mpNuclideLibrary->PartitionFunctionsWithoutSpinTerm(
      thermoState.T9(), mOpts.PartitionFunctionLog10Cap);

  if (mOpts.EnableScreening && mEOSEnabled
      && (thermoState.T9() >= mOpts.MinT9ForScreening)) {
    auto screen = mpScreening->Compute(thermoState, mCurrentY,
        *mpNuclideLibrary);

    for (auto& lib : mpReactionLibraries) {
        lib->SetTime(time);
        lib->CalculateRates(thermoState, partf, mOpts.RateExpArgumentCap,
            &mpNuclideLibrary->Zs(), &screen.Mus());
      }
  } else {
    for (auto& lib : mpReactionLibraries) {
      lib->SetTime(time);
      lib->CalculateRates(thermoState, partf, mOpts.RateExpArgumentCap,
          nullptr, nullptr);
    }
  }

  mProfiler.RateCalculation.Stop();
}

void ReactionNetwork::AddYdotContributions(const std::vector<double>& Y,
    std::vector<double> * const pYdot) const {
  mProfiler.YdotCalculation.Start();

  for (auto& lib : mpReactionLibraries)
    lib->AddYdotContributions(Y, pYdot);

  mProfiler.YdotCalculation.Stop();
}

void ReactionNetwork::AddJacobianContributions(const std::vector<double>& Y,
    Jacobian * const pJ) const {
  for (auto& lib : mpReactionLibraries)
    lib->AddJacobianContributions(Y, pJ);
}

ReactionNetwork::DtStr ReactionNetwork::FindMaxAbsXByYAndLimitingNuclideName(
    const std::vector<double>& x) const {
  std::vector<double> xByY(mCurrentY.size());
#if defined(__ICC) || defined(__INTEL_COMPILER)
  // if the Intel compiler vectorizes this loop, it can trigger a
  // floating point exception when it divides by zero (even though
  // it would not use the result due to the if statement, the
  // floating point exception still triggers)
#pragma novector
#endif // Intel
  for (unsigned int i = 0; i < xByY.size(); ++i) {
    if (mCurrentY[i] != 0.0) {
      xByY[i] = x[i] / mCurrentY[i];
    } else {
      if (x[i] == 0.0)
        xByY[i] = 1.0;
      else
        xByY[i] = std::numeric_limits<double>::max();
    }
  }

  double maxAbsXByY = 0.0;
  std::string limitingNuclideName = "";
  for (unsigned int i = 0; i < xByY.size(); ++i) {
    //printf("%u: yOld = %g, yDotOld = %g\n", i, yOld(i), yDotOld(i));
    if (mCurrentY[i] >= mOpts.SmallestYUsedForDtCalculation) {
      if (fabs(xByY[i]) > maxAbsXByY) {
        maxAbsXByY = fabs(xByY[i]);
        limitingNuclideName = mpNuclideLibrary->Names()[i];
      }
    }
  }

  return DtStr(maxAbsXByY, limitingNuclideName);
}

ReactionNetwork::DtStr ReactionNetwork::CalculateFirstDt(
    const double dtMax, const double dtInit) {
  double dt = std::min(dtMax, mOpts.MaxDt);
  std::string limitStr = "dtMax";

  double useDtInit = dtInit;
  if (dtInit <= 0.0)
    useDtInit = 1.0E4 * mOpts.MinDt;

  if (dtInit < dt) {
    dt = useDtInit;
    limitStr = "dtIni";
  }

  UpdateRates(mCurrentTime);
  std::vector<double> yDot(mCurrentY.size());

  AddYdotContributions(mCurrentY, &yDot);

  // find new time step
  DtStr maxYDotByYAndLimitingNuclideName =
      FindMaxAbsXByYAndLimitingNuclideName(yDot);
  double maxYDotByY = maxYDotByYAndLimitingNuclideName.Dt;

  if (maxYDotByY > 0.0) {
    double dtFromMaxYDotByY = 0.01 * mOpts.MaxYChangePerStep / maxYDotByY;
    if (dtFromMaxYDotByY < dt) {
      dt = dtFromMaxYDotByY;
      limitStr = maxYDotByYAndLimitingNuclideName.Str;
    }
  }

  return DtStr(dt, limitStr);
}

ReactionNetwork::DtStr ReactionNetwork::CalculateDt(const double dtMax) const {
  if (mPreviousY.size() == 0) {
    // this is the first step and we don't have previous Ys
    return DtStr(std::min(dtMax, mPreviousDtAndLimitString.Dt),
        mPreviousDtAndLimitString.Str);
  }

  std::vector<double> dY(mCurrentY.size());
  for (unsigned int i = 0; i < dY.size(); ++i) {
    dY[i] = mCurrentY[i] - mPreviousY[i];
  }

  DtStr maxDYByYAndLimitingNuclideName =
      FindMaxAbsXByYAndLimitingNuclideName(dY);
  double maxDYByY = maxDYByYAndLimitingNuclideName.Dt;
  std::string limitNuclideName = maxDYByYAndLimitingNuclideName.Str;

  std::string limitStr = "dtMax";
  double dt = std::min(std::min(dtMax, mOpts.MaxDt),
      mPreviousDtAndLimitString.Dt * mOpts.MaxDtChangeMultiplier);
  if (maxDYByY > 0.0) {
    double dtFromMaxYChange =
        mOpts.MaxYChangePerStep * mPreviousDtAndLimitString.Dt / maxDYByY;
    if (dtFromMaxYChange < dt) {
      dt = dtFromMaxYChange;
      limitStr = limitNuclideName;
    }
  }

  return DtStr(dt, limitStr);
}

void ReactionNetwork::UpdateRates(const double trialTime) {
  double rho = (*mpRhoVsTime)(trialTime);
  auto state = mCurrentThermodynamicState;

  // TODO reduce code duplication
  if (mEOSEnabled) {
    try {
      if (mOpts.IsSelfHeating) {
        // get temperature from entropy at beginning of the time step but density
        // at the end of the step
        state = mpEOS->FromSAndRhoWithTGuess(state.S(), rho, mCurrentY,
            *mpNuclideLibrary, state.T9(), (*mpNuDistVsTime)(trialTime));
      } else {
        state = mpEOS->FromTAndRho((*mpT9VsTime)(trialTime), rho, mCurrentY,
            *mpNuclideLibrary, (*mpNuDistVsTime)(trialTime));
      }
    } catch (std::exception& ex) {
      if (state.T9() >= mOpts.MinT9ForStrongReactions) {
        // this is a problem, we have to abort
        throw std::runtime_error("Exception in EOS in UpdateRates: "
            + std::string(ex.what()));
      } else {
        // we don't need the temperature anymore, just turn off EOS
        mEOSEnabled = false;
        mpOutput->Log("# Freezing temperature (EOS inversion failed)\n");
      }
    }
  }

  if (!mEOSEnabled) {
    mCurrentThermodynamicState.SetRho(rho);
    mCurrentThermodynamicState.SetNeutrinoDistribution(
        (*mpNuDistVsTime)(trialTime));
  }

  CalculateRates(trialTime, state);
}

HeatEntropy ReactionNetwork::CalculateHeatingRateAndEntropyChangeOrDerivative(
    const double time, const double dtOfLastStep, const double T9,
    const double rho, const std::vector<double>& currentY,
    const std::vector<double>& dYOrPreviousY,
    const bool calculateEntropyDerivative) const {
  HeatEntropy result;
  result.HeatingRate = 0.0;
  result.EntropyChangeOrDerivative = 0.0;

  if ((dtOfLastStep > 0.0) || calculateEntropyDerivative) {
    // calculate heating rate and change in entropy (or entropy derivative),
    // no matter whether the network is self-heating or not
    double constantLogContribution =
        1.5 * log(pow(rho, 2.0 / 3.0) * Constants::TwoPiHbar2C2NA23 / T9);

    double sumDYMassExcess = 0.0;
    double sumDYLogTerm = 0.0;
    double sumDYZ = 0.0;

    auto gs = mpNuclideLibrary->PartitionFunctionsWithoutSpinTerm(T9,
        mOpts.PartitionFunctionLog10Cap);

    for (unsigned int i = 0; i < currentY.size(); ++i) {
      double Y = currentY[i];

      if (Y < mOpts.SmallestYUsedForEntropyCalculation)
        continue;

      double dY = 0.0;
      if (calculateEntropyDerivative)
        dY = dYOrPreviousY[i]; // this is dY/dt
      else
        dY = Y - dYOrPreviousY[i]; // this is \Delta Y

      double G = mpNuclideLibrary->GroundStateSpinTerms()[i] * gs[i];
      double logTerm = log(Y / G)
          - 1.5 * log(mpNuclideLibrary->MassesInMeV()[i])
          + constantLogContribution;

      sumDYLogTerm += dY * logTerm;
      sumDYMassExcess += dY * mpNuclideLibrary->MassExcessesInMeV()[i];
      sumDYZ += dY * (double)mpNuclideLibrary->Zs()[i];
    }

    double TMeV = T9 * Constants::BoltzmannConstantInMeVPerGK;
    double gamma = Constants::ElectronMassInMeV / TMeV;
    // entropy change due to electrons
    double dSElectron = (mCurrentThermodynamicState.EtaElectron() + gamma)
        * sumDYZ;

    double dotEpsilonNu = 0.0; // in MeV per second
    for (auto& lib : mpReactionLibraries) {
      dotEpsilonNu += lib->GetEDotExternal(currentY);
    }

    double heatingRateInMeVPerSecPerBaryon = 0.0;

    if (calculateEntropyDerivative)
      heatingRateInMeVPerSecPerBaryon = -dotEpsilonNu - sumDYMassExcess;
    else
      heatingRateInMeVPerSecPerBaryon = -dotEpsilonNu
          - sumDYMassExcess / dtOfLastStep;

    if (mpExternalHeatingRate != nullptr)
      heatingRateInMeVPerSecPerBaryon += (*mpExternalHeatingRate)(time) /
          (Constants::ErgPerMeV * Constants::AvogadroConstantInPerGram);

    // convert heating rate from MeV / s / baryon to erg / s / g
    result.HeatingRate = heatingRateInMeVPerSecPerBaryon * Constants::ErgPerMeV
        * Constants::AvogadroConstantInPerGram;

    double currentTemperatureInMeV =
        T9 * Constants::BoltzmannConstantInMeVPerGK;
    if (calculateEntropyDerivative) {
      // return entropy derivative
      result.EntropyChangeOrDerivative = heatingRateInMeVPerSecPerBaryon
          / currentTemperatureInMeV - sumDYLogTerm - dSElectron;
    } else {
      // return entropy change
      result.EntropyChangeOrDerivative =
          heatingRateInMeVPerSecPerBaryon * dtOfLastStep
              / currentTemperatureInMeV - sumDYLogTerm - dSElectron;
    }
  }

  return result;
}

void ReactionNetwork::PostStepUpdate() {
  // calculate Ye and total mass
  double totalMass = 0.0;
  for (unsigned int i = 0; i < mCurrentY.size(); ++i)
    totalMass += mCurrentY[i] * (double)mpNuclideLibrary->As()[i];

  std::pair<double, std::string> dtStrPair;
  dtStrPair.first = mPreviousDtAndLimitString.Dt;
  dtStrPair.second = mPreviousDtAndLimitString.Str;
  mpOutput->AddEntry(mCurrentTime, dtStrPair, mCurrentNumberOfFailedTimeSteps,
      mCurrentNumberOfNRIterations, mCurrentThermodynamicState,
      mCurrentHeatingRate, mCurrentY, totalMass);

  // update state of reaction libraries
  bool stateChanged = false;

  // Output rates to screen
  for (auto& lib : mpReactionLibraries) {
    if (lib->GetDoOutput())
      lib->PrintRates(mpOutput.get());
  }

  for (auto& lib : mpReactionLibraries)
    stateChanged |= lib->StateChanged(mCurrentThermodynamicState,
        mpOutput.get());

  if (stateChanged) {
    mpOutput->Log("# Re-initializing Jacobian\n");

    // now the sparsity structure of the Jacobian has potentially changed,
    // we re-initialize the Jacobian with the new sparsity structure
    InitializeJacobian();
  }

  if (mEOSEnabled && mOpts.IsSelfHeating) {
    if (mCurrentThermodynamicState.Rho() < mOpts.MinDensityForEOS) {
      mEOSEnabled = false;
      mpOutput->Log("# Freezing temperature (density below threshold)\n");
    }
  }
}

double ReactionNetwork::TakeStep(const double dtInit,
    const bool calledRecursively) {
  std::vector<double> yNew(mCurrentY.size());
  double dt = std::max(dtInit, mOpts.MinDt);

  mCurrentNumberOfFailedTimeSteps = 0;

  while (!TryTakeStep(dt, &yNew)) {
    // reduce time step and try again
    dt /= mOpts.MaxDtChangeMultiplier;
    ++mCurrentNumberOfFailedTimeSteps;

    if (dt < mOpts.MinDt) {
      // try switching to NSE evolution
      if (SwitchToNSEEvolutionIfPossible(
          mOpts.MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE,
          "time step too small", &yNew)) {
        dt = mOpts.MinDt;
      } else {
        if (!calledRecursively) {
          RenormalizeMass();
          mpOutput->Log("# Renormalizing mass (time step too small)\n");
          return TakeStep(dtInit, true);
        } else {
          throw std::runtime_error("Newton-Raphson did not converge for "
              "smallest possible time step.");
        }
      }

      break;
    }
  }

  mPreviousY.swap(mCurrentY);
  mCurrentY.swap(yNew);

  // check whether this was a stuck step
  if ((mCurrentNumberOfFailedTimeSteps == 1)
      && (mpOutput->DtVsTime()[mpOutput->NumEntries() - 1] == dt)) {
    // this is a stuck step
    ++mCurrentNumberOfStuckSteps;

    if (mCurrentNumberOfStuckSteps > mOpts.NumStuckStepsForRenormalization) {
      mpOutput->Log("# Renormalizing mass (number of stuck steps for "
          "renormalization exceeded)\n");
      RenormalizeMass();
    }
  } else {
    mCurrentNumberOfStuckSteps = 0;
  }

  return dt;
}

bool ReactionNetwork::TryTakeStep(const double dt,
    std::vector<double> * const pYNew) {
  double tTrail = mCurrentTime + dt;
  UpdateRates(tTrail);

  // initial guess
  *pYNew = mCurrentY;
  std::vector<double> yDotNew(pYNew->size(), 0.0);

  double lastError = std::numeric_limits<double>::max();
  mCurrentNumberOfNRIterations = 0;
  while (true) {
    if (mCurrentNumberOfNRIterations > mOpts.MaxIterationsPerStep)
      return false;

    yDotNew = std::vector<double>(pYNew->size(), 0.0);
    AddYdotContributions(*pYNew, &yDotNew);

    std::vector<double> rhs(mCurrentY.size());
    for (unsigned int i = 0; i < rhs.size(); ++i)
      rhs[i] = yDotNew[i] - ((*pYNew)[i] - mCurrentY[i]) / dt;

    mProfiler.JacobianCalculation.Start();
    mpJacobian->Reset();

    AddJacobianContributions(*pYNew, mpJacobian.get());

    for (unsigned int i = 0; i < rhs.size(); ++i) {
      (*mpJacobian)(i, i) -= 1.0 / dt;
    }
    mProfiler.JacobianCalculation.Stop();

    mProfiler.MatrixInversion.Start();
    std::vector<double> deltaY = mpMatrixSolver->Solve(mpJacobian.get(), rhs);
    mProfiler.MatrixInversion.Stop();

    for (unsigned int i = 0; i < deltaY.size(); ++i) {
      (*pYNew)[i] -= deltaY[i];

      if ((*pYNew)[i] < 0.0) {
        (*pYNew)[i] = 0.0;
      }
    }

    ++mTotalNumberOfNRIterations;
    ++mCurrentNumberOfNRIterations;

    // check whether we have converged
    double error = 0.0;
    double deltaYByYError = 0.0;
    double massError = 0.0;

    if ((mOpts.ConvergenceCriterion ==
        NetworkConvergenceCriterion::DeltaYByY) ||
        (mOpts.ConvergenceCriterion ==
            NetworkConvergenceCriterion::BothDeltaYByYAndMass)) {
      for (unsigned int i = 0; i < deltaY.size(); ++i) {
        if ((*pYNew)[i] >= mOpts.SmallestYUsedForErrorCalculation) {
          double err = fabs(deltaY[i] / (*pYNew)[i]);
          deltaYByYError = std::max(deltaYByYError, err);
        }
      }

      error = deltaYByYError;
    }

    if ((mOpts.ConvergenceCriterion ==
        NetworkConvergenceCriterion::Mass) ||
        (mOpts.ConvergenceCriterion ==
            NetworkConvergenceCriterion::BothDeltaYByYAndMass)) {
      for (unsigned int i = 0; i < deltaY.size(); ++i) {
        massError += fabs((*pYNew)[i]) * (double)mpNuclideLibrary->As()[i];
      }

      massError = fabs(massError - 1.0);
      error = massError;
    }

    if (mOpts.ConvergenceCriterion ==
        NetworkConvergenceCriterion::BothDeltaYByYAndMass) {
      error = std::max(deltaYByYError, massError);
    }

    bool converged = false;
    switch (mOpts.ConvergenceCriterion) {
    case NetworkConvergenceCriterion::DeltaYByY:
      converged = (deltaYByYError < mOpts.DeltaYByYThreshold);
      break;
    case NetworkConvergenceCriterion::Mass:
      converged = (massError < mOpts.MassDeviationThreshold);
      break;
    case NetworkConvergenceCriterion::BothDeltaYByYAndMass:
      converged = ((deltaYByYError < mOpts.DeltaYByYThreshold)
          && (massError < mOpts.MassDeviationThreshold));
      break;
    }

    if (converged) {
      // Newton-Raphson succeeded, check temperature and entropy didn't change
      // too much
      double newTime = mCurrentTime + dt;

      auto heatEntropy = CalculateHeatingRateAndEntropyChangeOrDerivative(
          newTime, dt, mCurrentThermodynamicState.T9(),
          mCurrentThermodynamicState.Rho(), *pYNew, mCurrentY, false);
      double newS = mCurrentThermodynamicState.S()
          + heatEntropy.EntropyChangeOrDerivative;

      double rho = (*mpRhoVsTime)(newTime);
      auto nuDist = (*mpNuDistVsTime)(newTime);

      // TODO reduce code duplication
      if (mEOSEnabled) {
        ThermodynamicState newState = mCurrentThermodynamicState;
        try {
          if (mOpts.IsSelfHeating) {
            newState = mpEOS->FromSAndRhoWithTGuess(newS, rho, *pYNew,
                *mpNuclideLibrary, newState.T9(), nuDist);

            if (newState.T9() >= mOpts.MinT9ForT9SChangeLimit) {
              // impose limit on temperature and entropy change
              if ((fabs(newState.T9() - mCurrentThermodynamicState.T9())
                  / mCurrentThermodynamicState.T9()
                  > mOpts.MaxFractionalTemperatureAndEntropyChange)
                  || (fabs(newState.S() - mCurrentThermodynamicState.S())
                      / mCurrentThermodynamicState.S()
                      > mOpts.MaxFractionalTemperatureAndEntropyChange)) {
                return false;
              }
            }
          } else {
            newState = mpEOS->FromTAndRho((*mpT9VsTime)(newTime), rho, *pYNew,
                *mpNuclideLibrary, nuDist);
          }
        } catch (std::exception& ex) {
          if (newState.T9() >= mOpts.MinT9ForStrongReactions) {
            // this is a problem, we have to abort
            throw std::runtime_error("Exception in EOS in TryTakeStep: "
                + std::string(ex.what()));
          } else {
            // we don't need the temperature anymore, just turn off EOS
            mEOSEnabled = false;
            mpOutput->Log("# Freezing temperature (EOS inversion failed)\n");
          }
        }
        mCurrentThermodynamicState = newState;
      }

      if (!mEOSEnabled) {
        mCurrentThermodynamicState.SetRho(rho);
        mCurrentThermodynamicState.SetS(newS);
        mCurrentThermodynamicState.SetNeutrinoDistribution(nuDist);
        if (!mOpts.IsSelfHeating) {
          mCurrentThermodynamicState.SetT9((*mpT9VsTime)(newTime));
        }
      }

      mCurrentHeatingRate = heatEntropy.HeatingRate;

      return true;
    } else if (error > lastError) {
      // error is growing
      return false;
      // TODO make an option for this threshold
    } else if (fabs(error - lastError) / lastError < 0.1) {
      // error did not decrease fast enough, abort
      return false;
    }

    lastError = error;
  }

  // this should never happen
  return false;
}

void ReactionNetwork::RenormalizeMass() {
  double totalMass = 0.0;
  for (unsigned int i = 0; i < mCurrentY.size(); ++i)
    totalMass += mCurrentY[i] * (double)mpNuclideLibrary->As()[i];

  if (totalMass != 0.0)
    for (unsigned int i = 0; i < mCurrentY.size(); ++i)
      mCurrentY[i] /= totalMass;
}
