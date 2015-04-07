/// \file ReactionNetwork.hpp
/// \author jlippuner
/// \since Jun 30, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONNETWORK_HPP_
#define SKYNET_REACTIONNETWORK_HPP_

#include <memory>
#include <vector>

#include "EquationsOfState/EOS.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/Screening.hpp"

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

#include "Network/NetworkOptions.hpp"
#include "Network/NetworkOutput.hpp"
#include "Network/NSE.hpp"
#include "Network/NSEOptions.hpp"
#include "Network/TemperatureDensityHistory.hpp"

#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"

#include "ODEIntegrators/AdaptiveODEIntegrator.hpp"
#include "ODEIntegrators/RungeKuttaMethods.hpp"

#include "Reactions/ReactionLibraryBase.hpp"

#include "Utilities/FunctionVsTime.hpp"
#include "Utilities/Profiler.hpp"

// C++11 does not include a function std::make_unique analogous to
// std::make_shared. According to Herb Sutter, this is "partly an oversight, and
// it will almost certainly be added in the future" (see
// http://herbsutter.com/gotw/_102). So we use the code for make_unique provided
// by Herb Sutter.
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

namespace { // unnamed so these things can only be used in this file

const double NaN = std::numeric_limits<double>::quiet_NaN();

struct HeatEntropy {
  double HeatingRate;
  double EntropyChangeOrDerivative;
};

} // namespace [unnamed]

struct OperatorSplitStepResult {
  OperatorSplitStepResult(const double heating, const double entropyChange) :
    Heating(heating),
    EntropyChange(entropyChange) {}

  double Heating;
  double EntropyChange;
};

typedef std::vector<const ReactionLibraryBase *> ReactionLibs;

class ReactionNetwork {
public:
  ReactionNetwork(const NuclideLibrary& nuclideLibrary,
      const ReactionLibs& reactionLibraries, const EOS * const pEOS,
      const Screening * const pScreening, const NetworkOptions& options,
      const NSEOptions& NSEOpts = NSEOptions());

  void LoadNeutrinoHistory(const std::shared_ptr<NeutrinoHistory>& nuHist) {
    mpNuDistVsTime = nuHist;
  }

  void LoadExternalHeatingRate(
      const FunctionVsTime<double>& externalHeatingRate) {
    mpExternalHeatingRate = externalHeatingRate.MakeUniquePtr();
  }

  const NuclideLibrary& GetNuclideLibrary() const {
    return *mpNuclideLibrary;
  }

  const NetworkOutput& Output() const {
    return *mpOutput;
  }

  const NetworkOutput& Evolve(const std::vector<double>& initialY,
      const double startTime, const double endTime,
      const FunctionVsTime<double> * const pT9VsTime,
      const FunctionVsTime<double> * const pRhoVsTime,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingWithInitialTemperature(
      const std::vector<double>& initialY, const double startTime,
      const double endTime, const double initialT9,
      const FunctionVsTime<double> * const pRhoVsTime,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingWithInitialEntropy(
      const std::vector<double>& initialY, const double startTime,
      const double endTime, const double initialEntropy,
      const FunctionVsTime<double> * const pRhoVsTime,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& Evolve(const std::vector<double>& initialY,
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingWithInitialTemperature(
      const std::vector<double>& initialY, const double initialT9,
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingWithInitialEntropy(
      const std::vector<double>& initialY, const double initialEntropy,
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveFromNSE(const double startTime,
      const double endTime, const FunctionVsTime<double> * const pT9VsTime,
      const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingFromNSEWithInitialTemperature(
      const double startTime, const double endTime, const double initialT9,
      const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingFromNSEWithInitialEntropy(
      const double startTime, const double endTime, const double initialEntropy,
      const FunctionVsTime<double> * const pRhoVsTime, const double Ye,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveFromNSE(
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingFromNSEWithInitialTemperature(
      const double initialT9,
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  const NetworkOutput& EvolveSelfHeatingFromNSEWithInitialEntropy(
      const double initialEntropy,
      const TemperatureDensityHistory& temperatureDensityHistoryFile,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  void InitOperatorSplitEvolution(const std::vector<double>& initialY,
      const double startTime, const double initialT9, const double initialRho,
      const std::string& outputPrefix, const double initDt = 0.0,
      std::vector<std::string> massFractionOutputNuclideNames = { });

  OperatorSplitStepResult OperatorSplitTakeStep(const double T9,
      const double rho, const double endTime);

  void FinishOperatorSplitEvolution();

private:
  enum class InitMode {
    Temperature,
    Entropy,
    TemperatureVsTime
  };

  struct YeEntropy {
    double Ye;
    double Entropy;

    YeEntropy& operator+=(const YeEntropy& other);
    YeEntropy operator*(const double a) const;
  };

  struct DtStr {
    DtStr(const double dt, const std::string& str) :
        Dt(dt),
        Str(str) {}

    double Dt;
    std::string Str;
  };

  ReactionNetwork(const NuclideLibrary& nuclideLibrary,
      const ReactionLibs& reactionLibraries, std::unique_ptr<EOS>&& pEOS,
      std::unique_ptr<Screening>&& pScreening,
      const NetworkOptions& networkOptions, const NSEOptions& nseOpts);

  void InitializeJacobian();

  const NetworkOutput& DoEvolution(const InitMode initMode,
      const std::vector<double>& initialY,
      const double startTime, const double endTime, const double initDt,
      const double initialT9, const double initialEntropy,
      const FunctionVsTime<double> * const pT9VsTime,
      const FunctionVsTime<double> * const pRhoVsTime,
      const double Ye, const std::string& outputPrefix,
      std::vector<std::string> massFractionOutputNuclideNames = {});

  void SetupEvolution(const InitMode initMode,
      const std::vector<double>& initialY,
      const double startTime, const double endTime, const double initDt,
      const double initialT9, const double initialEntropy,
      const FunctionVsTime<double> * const pT9VsTime,
      const FunctionVsTime<double> * const pRhoVsTime,
      const double Ye, const std::string& outputPrefix,
      std::vector<std::string> massFractionOutputNuclideNames);

  void DoEvolutionSteps(const double endTime);

  void FinishEvolution();

  void SetupNSEEvolution();

  double CheckIfNSEEvolution();

  bool SwitchToNSEEvolutionIfPossible(const double fractionalJumpLimit,
      const std::string& reason, std::vector<double> * const pY);

  void EvolveNSE(const double dtMax);

  void CalculateRates(const double time, const ThermodynamicState thermoState);

  void AddYdotContributions(const std::vector<double>& Y,
      std::vector<double> * const pYdot) const;

  void AddJacobianContributions(const std::vector<double>& Y,
      Jacobian * const pJ) const;

  DtStr FindMaxAbsXByYAndLimitingNuclideName(
      const std::vector<double>& x) const;

  DtStr CalculateFirstDt(const double dtMax, const double dtInit);

  DtStr CalculateDt(const double dtMax) const;

  /// Calculate the reaction rates at the given trail time (end of the next time
  /// step).
  void UpdateRates(const double trialTime);

  /// If calculateEntropyDerivative is false, then pYdot must be nullptr and
  /// this function calculates the heating rate of the last time step using
  /// the changes in the abundances and it calculates the change in entropy
  /// from the changes in the abundances. If calculateEntropyDerivative is true,
  /// pYdot must not be nullptr and this function calculates the heating rate
  /// from the given pYdot and calculates the derivative of the entropy (at
  /// constant temperature and density) calculated from the given pYdot.
  /// The returned array has the heating rate in the first element and the
  /// change in entropy or entropy derivative in the second element.
  HeatEntropy CalculateHeatingRateAndEntropyChangeOrDerivative(
      const double time, const double dtOfLastStep, const double T9,
      const double rho, const std::vector<double>& currentY,
      const std::vector<double>& dYOrPreviousY,
      const bool calculateEntropyDerivative) const;

  /// This function is called after a successful time step. If the network is
  /// self-heating, this function calculates the increase in entropy and the
  /// new temperature. This function also stores various quantities for output.
  void PostStepUpdate();

  /// Takes a single time step with the given dtInit, which is the initial dt to
  /// use. The dt may be reduced to ensure convergence. Returns the actual dt
  /// that was used.
  double TakeStep(const double dtInit, const bool calledRecursively);

  /// Tries to take a single time step with final time mT + dt, writing the new
  /// abundances in the given vector. Returns true if the step succeeded.
  bool TryTakeStep(const double dt, std::vector<double> * const pYNew);

  /// Renormalize the mass to get total mass exactly equal to 1.
  void RenormalizeMass();

  // these objects store the nuclear data and reaction data
  std::unique_ptr<NuclideLibrary> mpNuclideLibrary;
  std::vector<std::unique_ptr<ReactionLibraryBase>> mpReactionLibraries;

  // the equation of state, Jacobian, matrix solver and the output
  std::unique_ptr<EOS> mpEOS;
  std::unique_ptr<Screening> mpScreening;
  std::unique_ptr<NSE> mpNSE;
  std::unique_ptr<Jacobian> mpJacobian;
  std::unique_ptr<MatrixSolver> mpMatrixSolver;
  std::unique_ptr<NetworkOutput> mpOutput;

  // keep track of what is enabled
  bool mEOSEnabled;

  // used for NSE evolution
  std::unique_ptr<RKF45<YeEntropy>> mpNseEvolutionRKF45;
  std::unique_ptr<AdaptiveODEIntegrator<std::array<YeEntropy, 2>>>
      mpNseEvolutionIntegrator;
  bool mDoNSEEvolution;
  AdaptiveODEIntegratorResult<std::array<YeEntropy, 2>> mNseStep;

  // these functions provide the temperature T9 in Giga Kelvin, density rho in
  // Grams per Centimeter^3
  std::unique_ptr<FunctionVsTime<double>> mpT9VsTime;
  std::unique_ptr<FunctionVsTime<double>> mpRhoVsTime;

  // possible neutrino history
  std::shared_ptr<NeutrinoHistory> mpNuDistVsTime;

  // possible external heating rate in ergs / s / g
  std::unique_ptr<FunctionVsTime<double>> mpExternalHeatingRate;

  // the network options
  const NetworkOptions mOpts;
  NSEOptions mNSEOpts;

  // current state
  double mCurrentTime; // in seconds
  ThermodynamicState mCurrentThermodynamicState;
  double mCurrentHeatingRate; // in erg / s / g
  std::vector<double> mCurrentY; // current abundances

  DtStr mPreviousDtAndLimitString;
  std::vector<double> mPreviousY; // abundances from previous step

  int mTotalNumberOfNRIterations;
  int mCurrentNumberOfNRIterations;
  int mCurrentNumberOfFailedTimeSteps;
  int mCurrentNumberOfStuckSteps;
  mutable Profiler mProfiler;
};

#endif // SKYNET_REACTIONNETWORK_HPP_
