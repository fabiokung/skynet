/// \file HelmholtzEOS.hpp
/// \author jlippuner
/// \since Aug 6, 2014
///
/// \brief
///
///

#ifndef SKYNET_EQUATIONSOFSTATE_HELMHOLTZEOS_HPP_
#define SKYNET_EQUATIONSOFSTATE_HELMHOLTZEOS_HPP_

#include <string>

#include "EquationsOfState/EOS.hpp"
#include "Network/NetworkOptions.hpp"

class HelmholtzEOS : public EOS {
public:
  HelmholtzEOS(const std::string& tablePath,
      const NetworkOptions& options = NetworkOptions());

  virtual std::unique_ptr<EOS> MakeUniquePtr() const {
    return std::unique_ptr<EOS>(new HelmholtzEOS(*this));
  }

  void SetOptions(const NetworkOptions& options) {
    mOpts = options;
  }

  virtual void PrintInfo(NetworkOutput * const pOutput) const;

  ThermodynamicState FromTAndRho(const double T9, const double rho,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

  ThermodynamicState FromSAndRho(const double entropy, const double rho,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

  ThermodynamicState FromSAndRhoWithTGuess(const double entropy, const double rho,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      const double T9Guess,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

  // TODO test this
  ThermodynamicState FromTAndS(const double T9, const double entropy,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

  ThermodynamicState FromTAndSWithRhoGuess(const double T9, const double entropy,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      const double rhoGuess,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

private:
  ThermodynamicState CalcEntropyAndDerivative(const double T9, const double rho,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
      std::shared_ptr<NeutrinoDistribution> nuDist
      = std::shared_ptr<NeutrinoDistribution>(
          new DummyNeutrinoDistribution())) const;

  // returns the specific electron and photon entropy from the FORTRAN code in
  // k_b per baryon
  ThermodynamicState GetElectronPhotonEntropyAndDerivative(const double T9,
      const double rho, const double Abar, const double Zbar) const;

  std::string mSource;
  NetworkOptions mOpts;
};

#endif // SKYNET_EQUATIONSOFSTATE_HELMHOLTZEOS_HPP_
