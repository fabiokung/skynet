/// \file NSEOptions.hpp
/// \author jlippuner
/// \since Jul 31, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_NSEOPTIONS_HPP_
#define SKYNET_NETWORK_NSEOPTIONS_HPP_

#include <array>

struct NSEOptions {
  NSEOptions():
    MaxIterations(1000),
    MaxIterationsWithGuess(100),
    MaxScreeningIterations(20),
    StepSizeLimit(1.0E-5),
    ConvergenceThreshold(1.0E-12),
    ScreeningConvergenceThreshold(1.0E-12),
    BasisZN0{{-1, 1}},
    BasisZN1{{0, 1}},
    DoScreening(true) { }

  int MaxIterations;
  int MaxIterationsWithGuess;
  int MaxScreeningIterations;

  double StepSizeLimit;
  double ConvergenceThreshold;
  double ScreeningConvergenceThreshold;

  std::array<int, 2> BasisZN0;
  std::array<int, 2> BasisZN1;

  bool DoScreening;

  // for SWIG
  const std::array<int, 2>& GetBasisZN0() const {
    return BasisZN0;
  }

  const std::array<int, 2>& GetBasisZN1() const {
    return BasisZN1;
  }

  void SetBasisZN0(const std::array<int, 2>& value) {
    BasisZN0 = value;
  }

  void SetBasisZN1(const std::array<int, 2>& value) {
    BasisZN1 = value;
  }

  void SetDoScreening(const bool doScreening) {
    DoScreening = doScreening;
  }
};

#endif // SKYNET_NETWORK_NSEOPTIONS_HPP_
