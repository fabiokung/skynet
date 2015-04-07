/// \file EOS.hpp
/// \author jlippuner
/// \since Aug 5, 2014
///
/// \brief
///
///

#ifndef SKYNET_EQUATIONSOFSTATE_EOS_HPP_
#define SKYNET_EQUATIONSOFSTATE_EOS_HPP_

#include <limits>
#include <memory>
#include <vector>

class NetworkOptions;

#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "NuclearData/NuclideLibrary.hpp"

class ThermodynamicState {
public:
  ThermodynamicState() :
      mT9(std::numeric_limits<double>::signaling_NaN()),
      mRho(std::numeric_limits<double>::signaling_NaN()),
      mS(std::numeric_limits<double>::signaling_NaN()),
      mU(std::numeric_limits<double>::signaling_NaN()),
      mEtaElectron(std::numeric_limits<double>::signaling_NaN()),
      mYe(std::numeric_limits<double>::signaling_NaN()),
      mdSdT9(std::numeric_limits<double>::signaling_NaN()) {
  }

  ThermodynamicState(const double T9, const double rho, const double S,
      const double U, const double etaElectron, const double Ye,
      const double dSdT9, std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) :
      mT9(T9),
      mRho(rho),
      mS(S),
      mU(U),
      mEtaElectron(etaElectron),
      mYe(Ye),
      mdSdT9(dSdT9),
      mNuDist(nuDist) {
    mNuDist->SetLocalT9(mT9);
  }

  static double CalcYe(const std::vector<double>& Y,
      const NuclideLibrary& nuclib) {
    double ye = 0.0;
    for (unsigned int i = 0; i < Y.size(); ++i)
      ye += Y[i] * (double)nuclib.Zs()[i];

    return ye;
  }

  static ThermodynamicState FromT9Rho(const double T9, const double rho,
      const std::vector<double>& Y, const NuclideLibrary& nuclib,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) {
    ThermodynamicState state;
    state.mT9 = T9;
    state.mRho = rho;
    state.mYe = ThermodynamicState::CalcYe(Y, nuclib);
    state.mNuDist = nuDist;
    state.mNuDist->SetLocalT9(T9);
    return state;
  }

  static ThermodynamicState FromT9RhoS(const double T9, const double rho,
      const double S, const std::vector<double>& Y,
      const NuclideLibrary& nuclib,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) {
    auto state = ThermodynamicState::FromT9Rho(T9, rho, Y, nuclib, nuDist);
    state.mS = S;

    return state;
  }

  double T9() const {
    return mT9;
  }

  void SetT9(const double T9) {
    mT9 = T9;
  }

  double Rho() const {
    return mRho;
  }

  void SetRho(const double rho) {
    mRho = rho;
  }

  double S() const {
    return mS;
  }

  void SetS(const double S) {
    mS = S;
  }

  double U() const {
    return mU;
  }

  void SetU(const double U) {
    mU = U;
  }

  double EtaElectron() const {
    return mEtaElectron;
  }

  void SetEtaElectron(const double etaElectron) {
    mEtaElectron = etaElectron;
  }

  double Ye() const {
    return mYe;
  }

  void SetYe(const double Ye) {
    mYe = Ye;
  }

  double dSdT9() const {
    return mdSdT9;
  }

  void SetDSdT9(const double dSdT9) {
    mdSdT9 = dSdT9;
  }

  std::shared_ptr<NeutrinoDistribution> GetNeutrinoDistribution() const {
    return mNuDist;
  }

  void SetNeutrinoDistribution(
      const std::shared_ptr<NeutrinoDistribution>& nuDist) {
    mNuDist = nuDist;
  }

private:
  double mT9;     // temperature in GK
  double mRho;    // density in gm / cm^3
  double mS;      // entropy in k_b / baryon
  double mU;      // internal energy in erg / gm
  double mEtaElectron; // Electron degeneracy parameter
  double mYe;     // electron fraction

  // derivatives
  double mdSdT9;  // derivative of entropy w.r.t. temperature
                  // in k_b / baryon / GK
  std::shared_ptr<NeutrinoDistribution> mNuDist;
};

class EOS {
public:
  virtual ~EOS() {}

  virtual std::unique_ptr<EOS> MakeUniquePtr() const =0;

  virtual void SetOptions(const NetworkOptions& options) =0;

  virtual void PrintInfo(NetworkOutput * const pOutput) const =0;

  virtual ThermodynamicState FromTAndRho(const double T9, const double rho,
      const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const =0;

  virtual ThermodynamicState FromSAndRho(const double entropy, const double rho,
      const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const =0;

  /// Same as above, but an initial guess for the temperature is provided.
  virtual ThermodynamicState FromSAndRhoWithTGuess(const double entropy,
      const double rho, const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary, const double T9Guess,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const =0;

  virtual ThermodynamicState FromTAndS(const double T9, const double entropy,
      const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const =0;

  /// Same as above, but an initial guess for the density is provided.
  virtual ThermodynamicState FromTAndSWithRhoGuess(const double T9,
      const double entropy, const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary, const double rhoGuess,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const =0;
};

#endif // SKYNET_EQUATIONSOFSTATE_EOS_HPP_
