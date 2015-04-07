/// \file HelmholtzEOS.cpp
/// \author jlippuner
/// \since Aug 6, 2014
///
/// \brief
///
///

#include "EquationsOfState/HelmholtzEOS.hpp"

#include <cmath>
#include <functional>

#include "Utilities/BisectionMethod.hpp"
#include "Utilities/CodeError.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/Exp10.hpp"
#include "Utilities/FlushStdout.hpp"

namespace { // unnamed so that these things can only be used in this file

extern "C" void read_table_(const char * const helmTablePath,
    const int * const str_length);

extern "C" void helm_eos_wrap_(
    const double * const T,               // temperature in K
    const double * const rho,             // density in g / cm^3
    const double * const Abar,            // average mass number A
    const double * const Zbar,            // average atomic number Z
    // output specific electron and photon entropy in erg / g / K
    double * const electronPhotonEntropy,
    // output specific internal electron and photon energy in erg / g
    double * const electronPhotonInternalEnergy,
    // Electron degeneracy parameter, dimensionless units
    double * const electronEta,
    // derivative of the specific electron and photon entropy in erg / g / K^2
    double * const electronPhotonEntropyDerivativeWithTemperature,
    bool * const success_flag             // true if succeeded
  );

} // namespace [unnamed]

HelmholtzEOS::HelmholtzEOS(const std::string& tablePath,
    const NetworkOptions& options) :
    mSource(tablePath),
    mOpts(options) {
  // the Helmholtz data table has to be read only once, so do this at the
  // beginning
  int len = tablePath.length();
  read_table_(tablePath.c_str(), &len);
}

void HelmholtzEOS::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# EOS: Default SkyNet EOS (modified Timmes)\n");
  pOutput->Log("#   Source: %s\n", mSource.c_str());
}

ThermodynamicState HelmholtzEOS::FromTAndRho(const double T9, const double rho,
    const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  return CalcEntropyAndDerivative(T9, rho, Y, nuclideLibrary, nuDist);
}

ThermodynamicState HelmholtzEOS::FromSAndRho(const double targetEntropy,
    const double rho, const std::vector<double>& Y,
    const NuclideLibrary& nuclideLibrary,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  return FromSAndRhoWithTGuess(targetEntropy, rho, Y, nuclideLibrary, 1.0,
      nuDist);
}

ThermodynamicState HelmholtzEOS::FromSAndRhoWithTGuess(
    const double targetEntropy, const double rho, const std::vector<double>& Y,
    const NuclideLibrary& nuclideLibrary, const double T9Guess,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  const double tolerance = 1.0E-16;

  // first do a bisection
  auto f = [&](const double log10T) {
    //printf("log10T = %f\n", log10T);
    auto res = this->FromTAndRho(exp10(log10T - 9.0), rho, Y, nuclideLibrary);
    return res.S() / targetEntropy - 1.0;
  };

  double log10Tres = BisectionMethod::FindRoot(f, log10(T9Guess) + 9.0,
      3.01, 12.99, 100, tolerance).Value;

  // tried to supplement this with a Newton-Raphson iteration, but that did not
  // improve things and made the algorithm much less robust

  //TODO add option to turn warning on or off
//  double error = f(log10Tres);
//  if (fabs(error) > tolerance) {
//    printf("# WARNING: Bisection only converged to an error %.4E,\n"
//        "# but requested tolerance was %.4E.\n", error, tolerance);
//  }

  return FromTAndRho(exp10(log10Tres - 9.0), rho, Y, nuclideLibrary, nuDist);
}

ThermodynamicState HelmholtzEOS::FromTAndS(const double T9,
    const double targetEntropy, const std::vector<double>& Y,
    const NuclideLibrary& nuclideLibrary,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  return FromTAndSWithRhoGuess(T9, targetEntropy, Y, nuclideLibrary, 1.0E6,
      nuDist);
}

ThermodynamicState HelmholtzEOS::FromTAndSWithRhoGuess(const double T9,
    const double targetEntropy, const std::vector<double>& Y,
    const NuclideLibrary& nuclideLibrary, const double rhoGuess,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  const double tolerance = 1.0E-16;

  // first do a bisection
  auto f = [&](const double log10Rho) {
    //printf("log10T = %f\n", log10T);
    auto res = this->FromTAndRho(T9, exp10(log10Rho - 11.0),
        Y, nuclideLibrary, nuDist);
    return 1.0 - res.S() / targetEntropy;
  };

  // shift density up by 11 orders of magnitude, since the argument to the
  // BisectionMethod cannot be negative
  double log10Rhores = BisectionMethod::FindRoot(f,
      log10(rhoGuess + 11.0), 2.0, 25.99, 100, tolerance).Value;

  return FromTAndRho(T9, exp10(log10Rhores - 11.0), Y, nuclideLibrary, nuDist);
}

ThermodynamicState HelmholtzEOS::CalcEntropyAndDerivative(const double T9,
    const double rho, const std::vector<double>& Y,
    const NuclideLibrary& nuclideLibrary,
    std::shared_ptr<NeutrinoDistribution> nuDist ) const {
  ASSERT((int )Y.size() == nuclideLibrary.NumNuclides(), "Number of abundances "
      "not equal to the number of nuclides.");

  // calculate the entropy contributions from each nuclear species and add them
  // up, also calculate Abar and Zbar at the same time
  double ionEntropy = 0.0;
  double ionInternalEnergy = 0.0;
  double ionEntropyDerivativeWithTemperature = 0.0;
  double sumY = 0.0;
  double Abar = 0.0;
  double Ye = 0.0;

  // 3/2 ln(k_B / (2 Pi hbar^2 c^2))
  double constEntropyContirbution =
      1.5 * log(Constants::BoltzmannConstantDivBy2PiHbar2C2InPerMeVPerGKPerCm2);

  double nb = Constants::AvogadroConstantInPerGram * rho;
  double kt = Constants::BoltzmannConstantsInErgPerK * T9 * 1.0E9;

  auto gs = nuclideLibrary.PartitionFunctionsWithoutSpinTerm(T9,
      mOpts.PartitionFunctionLog10Cap);
  auto dLogGByLogTs =
      nuclideLibrary.LogDerivativeOfPartitionFunctionsWithSpinTerm(T9,
          mOpts.PartitionFunctionLog10Cap);

  for (unsigned int i = 0; i < Y.size(); ++i) {
    double y = Y[i];
    if ((y < 0.0) || (y > 1.0)) {
      char buf[256];
      sprintf(buf, "Invalid Y = %.3E for %s in Helmholtz EOS", y,
          nuclideLibrary.Names()[i].c_str());
      throw std::invalid_argument(std::string(buf));
    }

    int A = nuclideLibrary.As()[i];
    int Z = nuclideLibrary.Zs()[i];

    sumY += y;
    Abar += y * (double)A;
    Ye += y * (double)Z;

    ionInternalEnergy += y * (kt * (1.5 + dLogGByLogTs[i])
        - nuclideLibrary.BindingEnergiesInMeV()[i] * Constants::ErgPerMeV);

    if (y > mOpts.SmallestYUsedForEntropyCalculation) {
      ionEntropy += y * (2.5
          + log(nuclideLibrary.GroundStateSpinTerms()[i] * gs[i] / (nb * y))
          + constEntropyContirbution
          + 1.5 * log(T9)
          + 1.5 * log(nuclideLibrary.MassesInMeV()[i])
          + dLogGByLogTs[i]);

      // ignoring second derivative of partition function, this needs to be
      // fixed at some point, but will probably require better partition
      // function tables
      ionEntropyDerivativeWithTemperature +=
          y / T9 * (1.5 + dLogGByLogTs[i] /* + second deriv of logG by logT */);
    }
  }

  Abar /= sumY;
  double Zbar = Ye / sumY;

  auto electronPhoton = GetElectronPhotonEntropyAndDerivative(T9, rho, Abar,
      Zbar);

  double npMassDiffInErg = (nuclideLibrary.NeutronMassInMeV()
      - nuclideLibrary.ProtonMassInMeV()) * Constants::ErgPerMeV;
  ionInternalEnergy -= Ye * npMassDiffInErg;

  ionInternalEnergy += Ye * Constants::ElectronMassInMeV * Constants::ErgPerMeV;

  // convert internal energy per baryon to per gram
  ionInternalEnergy *= Constants::AvogadroConstantInPerGram;

  double s = electronPhoton.S() + ionEntropy;
  double u = electronPhoton.U() + ionInternalEnergy;
  double dsdt = electronPhoton.dSdT9() + ionEntropyDerivativeWithTemperature;

  return ThermodynamicState(T9, rho, s, u, electronPhoton.EtaElectron(), Ye,
      dsdt, nuDist);
}

ThermodynamicState HelmholtzEOS::GetElectronPhotonEntropyAndDerivative(
    const double T9, const double rho, const double Abar,
    const double Zbar) const {
  double temperatureInKelvin = T9 * 1.0E9;
  //printf("T9 = %.10E, T = %.10E\n", T9, temperatureInKelvin);

  double ZbarSafe = Zbar + 1.0E-100;
  double s, u, etaElectron;
  double dsdt;
  bool successFlag = false;
  helm_eos_wrap_(&temperatureInKelvin, &rho, &Abar, &ZbarSafe, &s, &u,
      &etaElectron, &dsdt, &successFlag);

  // change units of entropy from erg / g / K to k_b per baryon
  double fac = Constants::BoltzmannConstantsInErgPerK
      * Constants::AvogadroConstantInPerGram;
  s /= fac;
  dsdt /= fac;

  // change units from K^{-1} to GK^{-1}
  dsdt *= 1.0E9;

  if (!successFlag) {
    FlushStdout::Flush();
    throw std::runtime_error("Helmholtz EOS failed");
  }

  return ThermodynamicState(T9, rho, s, u, etaElectron, Zbar / Abar, dsdt);
}
