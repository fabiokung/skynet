/// \file NetworkOptions.hpp
/// \author jlippuner
/// \since Jul 7, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_NETWORKOPTIONS_HPP_
#define SKYNET_NETWORK_NETWORKOPTIONS_HPP_

#include "Network/NetworkOutput.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name NetworkConvergenceCriterion in Python instead of
// the enum values being global constants

#ifdef SWIG
%rename(NetworkConvergenceCriterion) NetworkConvergenceCriterionStruct;
#endif // SWIG

struct NetworkConvergenceCriterionStruct {
  enum Value {
    /// The Newton-Raphson iteration is considered converged when DeltaY / Y is
    /// smaller than DeltaYByYThreshold.
    DeltaYByY,

    /// The Newton-Raphson iteration is considered converged when the difference
    /// between the total mass and 1 is smaller than MassDeviationThreshold.
    Mass,

    /// Both the DeltaYByY and Mass convergence criterion must be satisfied for
    /// the Newton-Raphson iteration to be considered converged.
    BothDeltaYByYAndMass
  };
};

typedef NetworkConvergenceCriterionStruct::Value NetworkConvergenceCriterion;

class NetworkOptions {
public:
  /// Create NetworkOptions class with default options.
  NetworkOptions();

  /// Write the options to the log file.
  void WriteToLog(NetworkOutput * const pOutput) const;

  /// Checks that all options make sense and throws an exception if not.
  void Check() const;

  /// Maximum number of Newton-Raphson iterations per time step.
  int MaxIterationsPerStep;

  /// Targeted maximum change in Y (abundances) per time step when
  /// calculating a time step size.
  double MaxYChangePerStep;

  /// Smallest abundance that is taken into account when calculating a time step
  /// size.
  double SmallestYUsedForDtCalculation;

  /// Convergence criterion that should be used.
  NetworkConvergenceCriterion ConvergenceCriterion;

  /// Minimum abundance to be considered when calculating the error in the
  /// Newton-Raphson iteration (only used for the DeltaYByY convergence
  /// criterion).
  double SmallestYUsedForErrorCalculation;

  /// With the DeltaYByY convergence criterion, the Newton-Raphson iteration is
  /// considered to have converged when the error is smaller than this
  /// threshold.
  double DeltaYByYThreshold;

  /// With the Mass convergence criterion, the Newton-Raphson iteration is
  /// considered to have converged when the difference between the total mass
  /// and 1 is less than this threshold.
  double MassDeviationThreshold;

  /// Maximum factor by which the time step can be multiplied or divided by.
  double MaxDtChangeMultiplier;

  /// Maximum dt is given by multiplying time scale over which density changes
  /// by this number.
  double DensityTimeScaleMultiplier;

  /// Smallest allowed time step size.
  double MinDt;

  /// Largest allowed time step size.
  double MaxDt;

  /// The maximum allowed fractional change in the temperature and entropy
  /// during a time step when the network is self-heating.
  double MaxFractionalTemperatureAndEntropyChange;

  /// If the temperature is below this limit, the maximum fractional change in
  /// the temperature and entropy is not enforced.
  double MinT9ForT9SChangeLimit;

  /// True of the network is self-heating.
  bool IsSelfHeating;

  /// True if electron screening is enabled.
  bool EnableScreening;

  /// Minimum temperature for screening to be used.
  double MinT9ForScreening;

  /// Minimum abundance to be used when calculating the change in entropy of a
  /// self-heating network. And also the minimum abundance to be used when
  /// calculating the entropy contribution of a particular nuclear species.
  double SmallestYUsedForEntropyCalculation;

  /// Cap for the argument that can go into an exp() when calculating a
  /// rate. If the argument is bigger than this, the rate is set to
  /// exp(RateExpArgumentCap).
  double RateExpArgumentCap;

  /// If the temperature falls below this value (in Giga Kelvin) all strong
  /// reactions are turned off.
  double MinT9ForStrongReactions;

  /// If the temperature is below this value (in Giga Kelvin), all neutrino
  /// reaction rates are set to 0.
  double MinT9ForNeutrinoReactions;

  /// If the density falls below this value (in g / cm^3), the temperature is
  /// frozen and no more EOS look-ups are performed.
  double MinDensityForEOS;

  /// True if regular output to stdout should be disabled
  bool DisableStdoutOutput;

  /// Cap for the log10 of a partition function value.
  double PartitionFunctionLog10Cap;

  /// NSE evolution can only be used if the temperature is above this threshold.
  double NSEEvolutionMinT9;

  /// Network evolution will not be used above this temperature, instead NSE evolution
  /// will be used
  double NetworkEvolutionMaxT9;

  /// If the density time scale becomes larger than this factor times the
  /// longest strong interaction time scale, NSE evolution is turned on.
  double NSEEvolutionStartTimeScaleFactor;

  /// The maximum allowed error in the ODEs for the NSE evolution. If the error
  /// is larger than this, the time step is reduced.
  double NSEEvolutionMaxError;

  /// The minimum error in the ODEs for the NSE evolution. If the error
  /// is smaller than this, the time step is increased.
  double NSEEvolutionMinError;

  /// The maximum allowed fractional change in the temperature and entropy
  /// that allows flashing to NSE if the time step gets smaller than MinDt.
  double MaxFractionalTemperatureAndEntropyChangeForEmergencyNSE;

  /// The mass is renormalized after this many stuck steps occurred. A stuck
  /// step is a time step that had exactly one failed attempt that resulted in
  /// the same dt as the previous step.
  int NumStuckStepsForRenormalization;
};

#endif // SKYNET_NETWORK_NETWORKOPTIONS_HPP_
