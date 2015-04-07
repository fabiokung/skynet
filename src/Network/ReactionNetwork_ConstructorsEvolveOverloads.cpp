/// \file ReactionNetwork_public.cpp
/// \author jlippuner
/// \since Aug 29, 2014
///
/// \brief
///
///

#include "Network/ReactionNetwork.hpp"

#include <limits>

#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

ReactionNetwork::ReactionNetwork(const NuclideLibrary& nuclideLibrary,
    const ReactionLibs& reactionLibraries, const EOS * const pEOS,
    const Screening * const pScreening, const NetworkOptions& options,
    const NSEOptions& nseOpts) :
    ReactionNetwork(nuclideLibrary, reactionLibraries,
        pEOS != nullptr ? pEOS->MakeUniquePtr() : nullptr,
        pScreening != nullptr ? pScreening->MakeUniquePtr() : nullptr,
        options, nseOpts) {
}

ReactionNetwork::ReactionNetwork(
    const NuclideLibrary& nuclideLibrary, const ReactionLibs& reactionLibraries,
    std::unique_ptr<EOS>&& pEOS, std::unique_ptr<Screening>&& pScreening,
    const NetworkOptions& options, const NSEOptions& nseOpts) :
    mpNuclideLibrary(nullptr),
    mpReactionLibraries(reactionLibraries.size()),
    mpEOS(std::move(pEOS)),
    mpScreening(std::move(pScreening)),
    mpJacobian(nullptr),
    mpOutput(nullptr),
    mEOSEnabled(true),
    mpNseEvolutionRKF45(nullptr),
    mpNseEvolutionIntegrator(nullptr),
    mDoNSEEvolution(true),
    mNseStep(),
    mpT9VsTime(nullptr),
    mpRhoVsTime(nullptr),
    mpNuDistVsTime(nullptr),
    mpExternalHeatingRate(nullptr),
    mOpts(options),
    mNSEOpts(nseOpts),
    mCurrentTime(NaN),
    mCurrentY(0),
    mPreviousDtAndLimitString(0.0, "none"),
    mPreviousY(0),
    mTotalNumberOfNRIterations(0),
    mCurrentNumberOfNRIterations(0),
    mCurrentNumberOfFailedTimeSteps(0) {
  // enable floating point exceptions
  FloatingPointExceptions::Enable();

  // check options
  mOpts.Check();

  if (mpEOS.get() == nullptr)
    throw std::invalid_argument("Constructed a reaction network with no EOS");

  // copy all reaction libraries
  for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i)
    mpReactionLibraries[i] = reactionLibraries[i]->MakeUniquePtr();

  // first collect all the nuclide names used in the reactions
  std::unordered_set<std::string> usedNuclideNames;
  for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i) {
    for (auto& reaction : mpReactionLibraries[i]->Reactions().ActiveData())
      reaction.AddUsedNuclideNames(&usedNuclideNames);
  }

  // restrict the nuclide library to only those nuclides that are used in the
  // reactions
  mpNuclideLibrary = make_unique<NuclideLibrary>(
      NuclideLibrary::CreateRestrictedLibrary(nuclideLibrary,
          std::vector<std::string>(usedNuclideNames.begin(),
              usedNuclideNames.end())));

  // since nuclides may have been thrown out in the above statement, we need to
  // re-index the nuclides in the reaction library to reflect the new nuclide
  // library
  for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i)
    mpReactionLibraries[i]->ReIndexNuclides(*mpNuclideLibrary);

  // if we have an EOS, set the current options on it
  if (mpEOS)
    mpEOS->SetOptions(mOpts);

  // set screening
  if (!mOpts.EnableScreening)
    mpScreening = nullptr;
  else if (mpScreening == nullptr)
    throw std::invalid_argument("Constructed a reaction network with "
        "screening enabled but no screening class was provided");

  // get maximum Z for which we have to compute screening corrections
  if (mOpts.EnableScreening) {
    int maxZ = 0;

    for (int i = 0; i < mpNuclideLibrary->NumNuclides(); ++i)
      maxZ = std::max(maxZ, mpNuclideLibrary->Zs()[i]);

    for (unsigned int i = 0; i < mpReactionLibraries.size(); ++i) {
      int thisMaxZ = mpReactionLibraries[i]->MaxScreeningZ(
          mpNuclideLibrary->Zs());
      maxZ = std::max(maxZ, thisMaxZ);
    }

    mpScreening->SetMaxZ(maxZ);
  }

  mNSEOpts.SetDoScreening(mOpts.EnableScreening);

  mpNSE = std::unique_ptr<NSE>(new NSE(*mpNuclideLibrary.get(), mpEOS.get(),
      mOpts.EnableScreening ? mpScreening.get() : nullptr, mNSEOpts));

  // create matrix solver
  mpMatrixSolver = MatrixSolver::Create();

  InitializeJacobian();

  // Set default neutrino history
  mpNuDistVsTime = std::unique_ptr<NeutrinoHistory>(new DummyNeutrinoHistory());
}

void ReactionNetwork::InitializeJacobian() {
  std::set<std::array<int, 2>> indicesOfNonZeroJacobianEntries;

  for (auto& lib : mpReactionLibraries)
    lib->AddIndicesOfNonZeroJacobianEntries(&indicesOfNonZeroJacobianEntries);

  mpJacobian = Jacobian::Create(mpNuclideLibrary->NumNuclides(),
      indicesOfNonZeroJacobianEntries, mpMatrixSolver.get());
}

const NetworkOutput& ReactionNetwork::Evolve(
    const std::vector<double>& initialY,
    const double startTime, const double endTime,
    const FunctionVsTime<double> * const pT9VsTime,
    const FunctionVsTime<double> * const pRhoVsTime,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::TemperatureVsTime, initialY, startTime, endTime,
      initDt, NaN, NaN, pT9VsTime, pRhoVsTime, NaN, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveSelfHeatingWithInitialTemperature(
    const std::vector<double>& initialY, const double startTime,
    const double endTime, const double initialT9,
    const FunctionVsTime<double> * const pRhoVsTime,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::Temperature, initialY, startTime, endTime,
      initDt, initialT9, NaN, nullptr, pRhoVsTime, NaN, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveSelfHeatingWithInitialEntropy(
    const std::vector<double>& initialY, const double startTime,
    const double endTime, const double initialEntropy,
    const FunctionVsTime<double> * const pRhoVsTime,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::Entropy, initialY, startTime, endTime,
      initDt, NaN, initialEntropy, nullptr, pRhoVsTime, NaN, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::Evolve(
    const std::vector<double>& initialY,
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return Evolve(initialY, temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(),
      &temperatureDensityHistoryFile.TemperatureVsTime(),
      &temperatureDensityHistoryFile.DensityVsTime(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveSelfHeatingWithInitialTemperature(
    const std::vector<double>& initialY, const double initialT9,
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return EvolveSelfHeatingWithInitialTemperature(initialY,
      temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(), initialT9,
      &temperatureDensityHistoryFile.DensityVsTime(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveSelfHeatingWithInitialEntropy(
    const std::vector<double>& initialY, const double initialEntropy,
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return EvolveSelfHeatingWithInitialEntropy(initialY,
      temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(), initialEntropy,
      &temperatureDensityHistoryFile.DensityVsTime(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveFromNSE(const double startTime,
    const double endTime, const FunctionVsTime<double> * const pT9VsTime,
    const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::TemperatureVsTime, { }, startTime, endTime,
      initDt, NaN, NaN, pT9VsTime, pRhoVsTime, Ye, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput&
ReactionNetwork::EvolveSelfHeatingFromNSEWithInitialTemperature(
    const double startTime, const double endTime, const double initialT9,
    const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::Temperature, { }, startTime, endTime,
      initDt, initialT9, NaN, nullptr, pRhoVsTime, Ye, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput&
ReactionNetwork::EvolveSelfHeatingFromNSEWithInitialEntropy(
    const double startTime, const double endTime, const double initialEntropy,
    const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return DoEvolution(InitMode::Entropy, { }, startTime, endTime,
      initDt, NaN, initialEntropy, nullptr, pRhoVsTime, Ye, outputPrefix,
      massFractionOutputNuclideNames);
}

const NetworkOutput& ReactionNetwork::EvolveFromNSE(
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return EvolveFromNSE(temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(),
      &temperatureDensityHistoryFile.TemperatureVsTime(),
      &temperatureDensityHistoryFile.DensityVsTime(),
      temperatureDensityHistoryFile.Ye(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}

const NetworkOutput&
ReactionNetwork::EvolveSelfHeatingFromNSEWithInitialTemperature(
    const double initialT9,
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return EvolveSelfHeatingFromNSEWithInitialTemperature(
      temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(),
      initialT9, &temperatureDensityHistoryFile.DensityVsTime(),
      temperatureDensityHistoryFile.Ye(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}

const NetworkOutput&
ReactionNetwork::EvolveSelfHeatingFromNSEWithInitialEntropy(
    const double initialEntropy,
    const TemperatureDensityHistory& temperatureDensityHistoryFile,
    const std::string& outputPrefix, const double initDt,
    std::vector<std::string> massFractionOutputNuclideNames) {
  return EvolveSelfHeatingFromNSEWithInitialEntropy(
      temperatureDensityHistoryFile.StartTime(),
      temperatureDensityHistoryFile.EndTime(),
      initialEntropy, &temperatureDensityHistoryFile.DensityVsTime(),
      temperatureDensityHistoryFile.Ye(), outputPrefix, initDt,
      massFractionOutputNuclideNames);
}
