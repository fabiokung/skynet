/// \file SkyNetScreening.cpp
/// \author jlippuner
/// \since Feb 22, 2017
///
/// \brief
///
///

#include "EquationsOfState/SkyNetScreening.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/FunctionIntegrator.hpp"

#include <string>

namespace { // unnamed so this can only be used in this file

inline double myF(const double x) {
  if (x > 700.0) return 0;
  if (x < -700.0) return 1.0;
  return 1.0 / (exp(x) + 1.0);
}

std::function<double(double)> MakeIntegrand(const double gamma,
    const double eta) {
  return [=] (const double x) {
    double arg = gamma * sqrt(1.0 + x * x);
    double fElectron = myF(arg - eta - gamma);
    double fPositron = myF(arg + eta + gamma);
    return x * x
        * (fElectron * (1.0 - fElectron) + fPositron * (1.0 - fPositron));
  };
}

} // namespace [unnamed]

void SkyNetScreening::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# Screening: Default SkyNet screening\n");
}

ScreeningState SkyNetScreening::Compute(const ThermodynamicState thermoState,
    const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary) const {
//  printf("SkyNetScreening::Compute:\n");
//  printf(" T9 = %.10e\n", thermoState.T9());
//  printf("rho = %.10e\n", thermoState.Rho());
//  printf("  s = %.10e\n", thermoState.S());
//  printf("  u = %.10e\n", thermoState.U());
//  printf("eta = %.10e\n", thermoState.EtaElectron());
//  printf(" ye = %.10e\n", thermoState.Ye());
//  printf(" dS = %.10e\n", thermoState.dSdT9());
//  printf("\n");

  const double b = 0.860;
  const double ln25 = log(25.0);

  double sumY = 0.0;
  double sumZY = 0.0;
  double sumZ2Y = 0.0;
  double sumZ3bm1Y = 0.0;

  for (size_t i = 0; i < Y.size(); ++i) {
    sumY += Y[i];

    double Z = (double)nuclideLibrary.Zs()[i];

    sumZY += Z * Y[i];
    sumZ2Y += Z * Z * Y[i];
    sumZ3bm1Y += pow(Z, 3 * b - 1) * Y[i];
  }

  double Zbar = sumZY / sumY;

  double nb = thermoState.Rho() * Constants::AvogadroConstantInPerGram;
  double TMeV = thermoState.T9() * Constants::BoltzmannConstantInMeVPerGK;

  double lambda0 = 1.9370131349470739099E-19 * sqrt(nb * sumY) / pow(TMeV, 1.5);

  FunctionIntegrator integrator;
  double gamma = Constants::ElectronMassInMeV / TMeV;
  auto integrand = MakeIntegrand(gamma, thermoState.EtaElectron());
  double integral = integrator.Integrate(integrand, 0.0,
      std::numeric_limits<double>::infinity());
  double me3 = Constants::ElectronMassInMeV * Constants::ElectronMassInMeV
      * Constants::ElectronMassInMeV;
  double zeta = sqrt(
      (sumZ2Y + 1.3186846724626433276E31 * me3 / nb * integral) / sumY);

  double etab = sumZ3bm1Y / sumY
      / (pow(zeta, 3 * b - 2) * pow(Zbar, 2 - 2 * b));

  std::vector<double> mus(mMaxZ + 1);
  mus[0] = 0.0;

  // the different mu's here are actually \beta * \mu(Z)

  for (int z = 1; z <= mMaxZ; ++z) {
    double Zd = (double)z;

    double p = (zeta + Zd) * zeta * zeta * lambda0;
    double twoLnP = 2.0 * log(p);

    double fw = 0.5 * (tanh(-twoLnP - ln25) + 1.0);
    double muWeak = -0.5 * Zd * Zd * zeta * lambda0;

    double fs = 0.5 * (tanh(twoLnP - ln25) + 1.0);
    double muStrong = -pow(lambda0, 2.0 / 3.0)
        * (0.6240 * pow(Zd, 5.0 / 3.0) * pow(Zbar, 1.0 / 3.0)
            + 0.1971 * pow(Zd, 4.0 / 3.0) * pow(Zbar, 2.0 / 3.0)
            - 0.0374 * Zd * Zbar) + 9.0 / 16.0 * Zd / Zbar
        - 0.4600 * pow(Zd / Zbar, 2.0 / 3.0);

    double fi = 1.0 - fw - fs;
    double muIntermediate = -0.380 * pow(lambda0, b) * etab * pow(Zd, b + 1);

    mus[z] = fw * muWeak + fi * muIntermediate + fs * muStrong;

//    printf("%3i  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e  %10.4e\n",
//        z, mus[z], muStrong, muWeak, muIntermediate, fs, fw, fi);
  }

  return ScreeningState(mus);
}



