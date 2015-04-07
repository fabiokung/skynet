/// \file NetworkOptions.cpp
/// \author jlippuner
/// \since Aug 25, 2014
///
/// \brief
///
///

#include "Network/NetworkOptions.hpp"

NetworkOptions::NetworkOptions() :
    MaxIterationsPerStep(10),
    MaxYChangePerStep(0.1),
    SmallestYUsedForDtCalculation(1.0E-6),
    ConvergenceCriterion(NetworkConvergenceCriterion::Mass),
    SmallestYUsedForErrorCalculation(1.0E-20),
    DeltaYByYThreshold(1.0E-10),
    MassDeviationThreshold(1.0E-10),
    MaxDtChangeMultiplier(2.0),
    DensityTimeScaleMultiplier(0.01),
    MinDt(1.0E-16),
    MaxDt(1.0E16),
    MaxFractionalTemperatureAndEntropyChange(0.01),
    MinT9ForT9SChangeLimit(1.0),
    IsSelfHeating(false),
    EnableScreening(false),
    MinT9ForScreening(1.0E-4),
    SmallestYUsedForEntropyCalculation(1.0E-20),
    RateExpArgumentCap(100.0),
    MinT9ForStrongReactions(0.01),
    MinT9ForNeutrinoReactions(0.001),
    MinDensityForEOS(1.0E-8),
    DisableStdoutOutput(false),
    PartitionFunctionLog10Cap(100.0),
    NSEEvolutionMinT9(7.0),
    NetworkEvolutionMaxT9(1.0E16),
    NSEEvolutionStartTimeScaleFactor(5.0),
    NSEEvolutionMaxError(1.0E-8),
    NSEEvolutionMinError(1.0E-10),
    MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE(0.5),
    NumStuckStepsForRenormalization(25) {}

void NetworkOptions::WriteToLog(NetworkOutput * const pOutput) const {
  pOutput->Log("#                 max iterations per step = %i\n",
      MaxIterationsPerStep);
  pOutput->Log("#                   max Y change per step = %12.6E\n",
      MaxYChangePerStep);
  pOutput->Log("#      smallest Y used for dt calculation = %12.6E\n",
      SmallestYUsedForDtCalculation);
  pOutput->Log("#                   convergence criterion = %s\n",
      (ConvergenceCriterion == NetworkConvergenceCriterion::DeltaYByY) ?
          "delta Y / Y small" :
      (ConvergenceCriterion == NetworkConvergenceCriterion::Mass) ?
          "mass conservation" :
      (ConvergenceCriterion ==
          NetworkConvergenceCriterion::BothDeltaYByYAndMass) ?
          "both delta Y / Y small and mass conservation" :
          throw std::invalid_argument("unknown convergence criterion"));
  pOutput->Log("#   smallest Y used for error calculation = %12.6E\n",
      SmallestYUsedForErrorCalculation);
  pOutput->Log("#                   delta Y / Y threshold = %12.6E\n",
      DeltaYByYThreshold);
  pOutput->Log("#                mass deviation threshold = %12.6E\n",
      MassDeviationThreshold);
  pOutput->Log("#                max dt change multiplier = %12.6E\n",
      MaxDtChangeMultiplier);
  pOutput->Log("#           density time scale multiplier = %12.6E\n",
      DensityTimeScaleMultiplier);
  pOutput->Log("#                                  min dt = %12.6E\n",
      MinDt);
  pOutput->Log("#                                  max dt = %12.6E\n",
      MaxDt);
  pOutput->Log("# max fractional change in T9 and entropy = %12.6E\n",
      MaxFractionalTemperatureAndEntropyChange);
  pOutput->Log("#        min T9 for T9 and S change limit = %12.6E\n",
      MinT9ForT9SChangeLimit);
  pOutput->Log("#                         is self-heating = %s\n",
      IsSelfHeating ? "true" : "false");
  pOutput->Log("#                               screening = %s\n",
      EnableScreening ? "true" : "false");
  pOutput->Log("#                    min T9 for screening = %12.6E\n",
        MinT9ForScreening);
  pOutput->Log("# smallest Y used for entropy calculation = %12.6E\n",
      SmallestYUsedForEntropyCalculation);
  pOutput->Log("#                   rate exp argument cap = %12.6E\n",
      RateExpArgumentCap);
  pOutput->Log("#             min T9 for strong reactions = %12.6E\n",
      MinT9ForStrongReactions);
  pOutput->Log("#           min T9 for neutrino reactions = %12.6E\n",
      MinT9ForNeutrinoReactions);
  pOutput->Log("#                     min density for EOS = %12.6E\n",
      MinDensityForEOS);
  pOutput->Log("#               stdout output is disabled = %s\n",
      DisableStdoutOutput ? "true" : "false");
  pOutput->Log("#            partition function log10 cap = %12.6E\n",
      PartitionFunctionLog10Cap);
  pOutput->Log("#                min T9 for NSE evolution = %12.6E\n",
      NSEEvolutionMinT9);
  pOutput->Log("#   NSE evolution start time scale factor = %12.6E\n",
      NSEEvolutionStartTimeScaleFactor);
  pOutput->Log("#             max error for NSE evolution = %12.6E\n",
      NSEEvolutionMaxError);
  pOutput->Log("#             min error for NSE evolution = %12.6E\n",
      NSEEvolutionMinError);
  pOutput->Log("# max frac chng in T9, s for small dt NSE = %12.6E\n",
      MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE);
  pOutput->Log("#  num of stuck steps for renormalization = %i\n",
      NumStuckStepsForRenormalization);
}

void NetworkOptions::Check() const {
  if (MaxIterationsPerStep < 1)
    throw std::invalid_argument("NetworkOptions.MaxIterationsPerStep must be "
        "at least 1.");

  if (MaxYChangePerStep <= 0.0)
    throw std::invalid_argument("NetworkOptions.MaxYChangePerStep must be "
        "positive.");

  if (SmallestYUsedForDtCalculation <= 0.0)
    throw std::invalid_argument("NetworkOptions.SmallestYUsedForDtCalculation "
        "must be positive.");

  if (SmallestYUsedForErrorCalculation <= 0.0)
    throw std::invalid_argument(
        "NetworkOptions.SmallestYUsedForErrorCalculation "
            "must be positive.");

  if (DeltaYByYThreshold <= 0.0)
    throw std::invalid_argument("NetworkOptions.DeltaYByYThreshold must be "
        "positive.");

  if (MassDeviationThreshold <= 0.0)
    throw std::invalid_argument("NetworkOptions.MassDeviationThreshold must be "
        "positive.");

  if (MaxDtChangeMultiplier <= 1.0)
    throw std::invalid_argument("NetworkOptions.MaxDtChangeMultiplier must be "
        "strictly greater than 1.");

  if (DensityTimeScaleMultiplier <= 0.0)
    throw std::invalid_argument("NetworkOptions.DensityTimeScaleMultiplier "
        "must be positive.");

  if (MinDt <= 0.0)
    throw std::invalid_argument("NetworkOptions.MinDt must be positive.");

  if (MaxDt <= MinDt)
    throw std::invalid_argument("NetworkOptions.MaxDt must be strictly greater "
        "than NetworkOptions.MinDt.");

  if (MaxFractionalTemperatureAndEntropyChange <= 0.0)
    throw std::invalid_argument("NetworkOptions."
        "MaxFractionalTemperatureAndEntropyChange must be positive.");

  if (MinT9ForT9SChangeLimit <= 0.0)
    throw std::invalid_argument("NetworkOptions.MinT9ForT9SChangeLimit must be "
        "positive.");

  if (SmallestYUsedForEntropyCalculation <= 0.0)
    throw std::invalid_argument("NetworkOptions."
        "SmallestYUsedForEntropyCalculation must be positive.");

  if (RateExpArgumentCap <= 0.0)
    throw std::invalid_argument("NetworkOptions.RateExpArgumentCap must be "
        "positive.");

  if (PartitionFunctionLog10Cap <= 0.0)
    throw std::invalid_argument("NetworkOptions.PartitionFunctionLog10Cap must "
        "be positive.");

  if (NSEEvolutionStartTimeScaleFactor <= 0.0)
    throw std::invalid_argument("NetworkOptions."
        "NSEEvolutionStartTimeScaleFactor must be positive.");

  if (NSEEvolutionMinError <= 0.0)
    throw std::invalid_argument("NetworkOptions.NSEEvolutionMinError must be "
        "positive.");

  if (NSEEvolutionMaxError <= NSEEvolutionMinError)
    throw std::invalid_argument("NetworkOptions.NSEEvolutionMaxError must be "
        "strictly greater than NSEEvolutionMinError.");

  if (MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE <= 0.0)
    throw std::invalid_argument("NetworkOptions."
        "MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE must be "
        "positive.");
}
