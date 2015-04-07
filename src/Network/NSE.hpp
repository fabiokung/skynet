/// \file NSE.hpp
/// \author jlippuner
/// \since Jul 30, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_NSE_HPP_
#define SKYNET_NETWORK_NSE_HPP_

#include <array>
#include <vector>

#include "EquationsOfState/EOS.hpp"
#include "EquationsOfState/Screening.hpp"
#include "Network/NSEOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"

namespace { // unnamed so this enum can only be used in this file

enum class NSEMode {
  FixedT9AndRho,
  FixedT9,
  FixedRho
};

class NSENotConverged : public std::exception {};

} // namespace [unnamed]

struct MuLM {
  double muL;
  double muM;

  MuLM() :
      muL(std::numeric_limits<double>::quiet_NaN()),
      muM(std::numeric_limits<double>::quiet_NaN()) {}
};

class NSE;

class NSEResult {
  friend NSE;
public:
  NSEResult(const std::vector<double>& y, const double T9, const double rho,
      const int numIterations, const double error, const MuLM muLM) :
      mY(y),
      mT9(T9),
      mRho(rho),
      mNumIterations(numIterations),
      mError(error),
      mMuLM(muLM) {}

  const std::vector<double>& Y() const {
    return mY;
  }

  double T9() const {
    return mT9;
  }

  double Rho() const {
    return mRho;
  }

  int NumIterations() const {
    return mNumIterations;
  }

  double Error() const {
    return mError;
  }

  MuLM GetMuLM() const {
    return mMuLM;
  }

private:
  std::vector<double> mY;
  double mT9;
  double mRho;
  int mNumIterations;
  double mError;
  MuLM mMuLM;
};

class NSE {
public:
  NSE(const NuclideLibrary& nuclib, const EOS * const pEOS,
      const Screening * const pScreening,
      const NSEOptions& opts = NSEOptions());

  NSEResult CalcFromTemperatureAndDensity(const double T9,
      const double rho, const double Ye);

  NSEResult CalcFromTemperatureAndDensityWithGuess(const double T9,
      const double rho, const double Ye, const MuLM muLMGuess);

  // TODO test this
  NSEResult CalcFromTemperatureAndEntropy(const double T9,
      const double entropy, const double Ye);

  // TODO test this
  NSEResult CalcFromTemperatureAndEntropyWithGuess(const double T9,
      const double entropy, const double Ye, const double rhoGuess,
      const MuLM muLMGuess);

  NSEResult CalcFromEntropyAndDensity(const double entropy,
      const double rho, const double Ye);

  NSEResult CalcFromEntropyAndDensityWithTemperatureGuess(
      const double entropy, const double rho, const double Ye,
      const double T9Guess);

  NSEResult CalcFromEntropyAndDensityWithGuess(const double entropy,
      const double rho, const double Ye, const double T9Guess,
      const MuLM muLMGuess);

  NSEResult CalcFromInternalEnergyAndDensity(const double internalEnergy,
      const double rho, const double Ye);

  NSEResult CalcFromInternalEnergyAndDensityWithTemperatureGuess(
      const double internalEnergy, const double rho, const double Ye,
      const double T9Guess);

  const NuclideLibrary& GetNuclideLibrary() const {
    return mNuclib;
  }

  const NSEOptions& GetNSEOptions() const {
    return mOpts;
  }

private:
  MuLM CalcInitialGuess(const double T9, const double rho,
      const double Ye);

  void CalcY(const double T9, const double etaL, const double etaM,
      const double rho, const int n, const std::vector<double>& alpha,
      const std::vector<double>& beta, std::vector<double> * const pY,
      const std::vector<double> * const pScreenCorr);

  double CalcF(const NSEMode mode, const double T9, const double rho,
      const double Ye, const double entropy, const int n,
      const std::vector<double>& Y, std::vector<double> * const pFs);

  void CalcJacobian(const NSEMode mode, const double T9,
      const double etaL, const double etaM, const double rho,
      const double entropy, const int n, const std::vector<double>& alpha,
      const std::vector<double>& beta, const std::vector<double>& Y,
      std::vector<std::vector<double>> * const pJ);

  NSEResult DoNRIteration(const NSEMode mode, const double inputT9,
      const double inputRho, const double Ye, const double entropy, const int n,
      const std::vector<double>& alpha, const std::vector<double>& beta,
      MuLM initialMuLM, const int maxIterations,
      const std::vector<double> * const pScreenCorr);

  NSEResult CalcNSE(const NSEMode mode, const double inputT9,
      const double inputRho, const double Ye, const double entropy,
      const MuLM muLMGuess);

  void CheckInput(const double T9, const double rho, const double Ye);

  const NuclideLibrary& mNuclib;
  std::unique_ptr<EOS> mpEOS;
  std::unique_ptr<Screening> mpScreening;
  NSEOptions mOpts;
};

#endif // SKYNET_NETWORK_NSE_HPP_
