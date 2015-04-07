/// \file NSE.cpp
/// \author jlippuner
/// \since Jul 30, 2014
///
/// \brief
///
///

#include "Network/NSE.hpp"

#include <cmath>
#include <stdexcept>

#include "Utilities/BisectionMethod.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/Exp10.hpp"
#include "Utilities/FloatingPointComparison.hpp"

namespace { // unnamed so this function can only be used in this file

double GetTempDensityFactor(const double TMeV, const double rho) {
  double nb = Constants::AvogadroConstantInPerGram * rho;
  //printf("TMeV = %.10E, nb = %.10E\n", TMeV, nb);
  return pow(TMeV / (2.0 * Constants::Pi
      * Constants::ReducedPlanckConstantInMeVSec
      * Constants::ReducedPlanckConstantInMeVSec), 3.0 / 2.0)
      / (nb * Constants::SpeedOfLightInCmPerSec
          * Constants::SpeedOfLightInCmPerSec
          * Constants::SpeedOfLightInCmPerSec);
}

inline double Det(
    const double a, const double b, const double c,
    const double d, const double e, const double f,
    const double g, const double h, const double i) {
  return a * (e * i - h * f) - b * (d * i - g * f) + c * (d * h - g * e);
}

struct T9Rho {
  double T9;
  double Rho;
};

T9Rho GetT9Rho(const NSEMode mode, const double inputT9, const double inputRho,
    const double T9OrRho) {
  T9Rho t9Rho;
  if (mode == NSEMode::FixedT9AndRho) {
    t9Rho.T9 = inputT9;
    t9Rho.Rho = inputRho;
  } else if (mode == NSEMode::FixedT9) {
    t9Rho.T9 = inputT9;
    t9Rho.Rho = T9OrRho;
  } else if (mode == NSEMode::FixedRho) {
    t9Rho.T9 = T9OrRho;
    t9Rho.Rho = inputRho;
  } else {
    throw std::runtime_error("Invalid NSEMode");
  }

  return t9Rho;
}

} // namespace [unnamed]

NSE::NSE(const NuclideLibrary& nuclib, const EOS * const pEOS,
    const Screening * const pScreening, const NSEOptions& opts) :
    mNuclib(nuclib),
    mpEOS(pEOS != nullptr ? pEOS->MakeUniquePtr() : nullptr),
    mpScreening(pScreening != nullptr ? pScreening->MakeUniquePtr() : nullptr),
    mOpts(opts) {
  if (mpEOS.get() == nullptr)
    throw std::invalid_argument("Constructed an NSE class with no EOS");

  if (mOpts.DoScreening && (mpScreening.get() == nullptr))
    throw std::invalid_argument("Screening enabled in NSE, but no screening "
        "class given");

  if (mOpts.DoScreening) {
    int maxZ = 0;
    for (int i = 0; i < nuclib.NumNuclides(); ++i)
      maxZ = std::max(maxZ, nuclib.Zs()[i]);

    if (mpScreening->MaxZ() < maxZ)
      mpScreening->SetMaxZ(maxZ);
  }
}

NSEResult NSE::CalcFromTemperatureAndDensity(const double T9,
    const double rho, const double Ye) {
  return CalcFromTemperatureAndDensityWithGuess(T9, rho, Ye, MuLM());
}

NSEResult NSE::CalcFromTemperatureAndDensityWithGuess(const double T9,
    const double rho, const double Ye,  const MuLM muLMGuess) {
  return CalcNSE(NSEMode::FixedT9AndRho, T9, rho, Ye,
      std::numeric_limits<double>::quiet_NaN(), muLMGuess);
}

NSEResult NSE::CalcFromTemperatureAndEntropy(const double T9,
    const double entropy, const double Ye) {
  double rhoGuess = 1.0E6;
  return CalcFromTemperatureAndEntropyWithGuess(T9, entropy, Ye, rhoGuess,
      MuLM());
}

NSEResult NSE::CalcFromTemperatureAndEntropyWithGuess(const double T9,
    const double entropy, const double Ye, const double rhoGuess,
    const MuLM /*muLMGuess*/) {
  // we only have the temperature and entropy, so we have to do a bisection to
  // find the density that gives this entropy for the NSE abundances obtained
  // with that density
  auto f = [&] (const double log10Rho) {
    double rho = exp10(log10Rho);
    auto res = CalcFromTemperatureAndDensity(T9, rho, Ye);
    auto eosRes = mpEOS->FromTAndRho(T9, rho, res.Y(), mNuclib);
    return 1.0 - eosRes.S() / entropy;
  };

  BisectionReturnValue log10RhoRoot(-1.0, false);
  try {
    log10RhoRoot = BisectionMethod::FindRoot(f, log10(rhoGuess), 1.0, 12.4);
  } catch (std::exception& ex) {
    throw std::runtime_error("Bisection failed in "
        "NSE::CalcFromTemperatureAndEntropyWithGuess with T9 = "
        + std::to_string(T9) + ", entropy = " + std::to_string(entropy)
        + ", Ye = " + std::to_string(Ye) + "\nCause: " + ex.what());
  }
  auto res = CalcFromTemperatureAndDensity(T9, exp10(log10RhoRoot.Value), Ye);

  return NSEResult(res.Y(), T9, exp10(log10RhoRoot.Value), res.NumIterations(),
      res.Error(), res.GetMuLM());
}

NSEResult NSE::CalcFromEntropyAndDensity(const double entropy,
    const double rho, const double Ye) {
  return CalcFromEntropyAndDensityWithTemperatureGuess(entropy, rho, Ye, 10.0);
}

NSEResult NSE::CalcFromEntropyAndDensityWithTemperatureGuess(
    const double entropy, const double rho, const double Ye,
    const double T9Guess) {
  auto f = [&] (const double log10T9) {
    double T9 = exp10(log10T9);
    auto res = CalcFromTemperatureAndDensity(T9, rho, Ye);
    auto eosRes = mpEOS->FromTAndRho(T9, rho, res.Y(), mNuclib);
    return eosRes.S() / entropy - 1.0;
  };

  // first do a quick bisection to get a good initial guess for NR
  BisectionReturnValue log10T9Root(-1.0, false);
  try {
    log10T9Root = BisectionMethod::FindRoot(f, log10(T9Guess), log10(1.001),
        log10(199.99), 100, 1.0E-2);
  } catch (std::exception& ex) {
    throw std::runtime_error("Bisection failed in "
        "NSE::CalcFromEntropyAndDensityWithTemperatureGuess with entropy = "
        + std::to_string(entropy) + ", rho = " + std::to_string(rho)
        + ", Ye = " + std::to_string(Ye) + "\nCause: " + ex.what());
  }

  auto res = CalcFromTemperatureAndDensity(exp10(log10T9Root.Value), rho, Ye);

  if (!log10T9Root.IsBestPossibleValue) {
    // we have a chance of getting a better result
    try {
      res = CalcFromEntropyAndDensityWithGuess(entropy, rho, Ye, res.T9(),
          res.GetMuLM());
    } catch (...) {
      // the NR iteration failed, do the best bisection we can
      log10T9Root = BisectionMethod::FindRoot(f, log10T9Root.Value,
          log10(1.001), log10(199.99));
      res = CalcFromTemperatureAndDensity(exp10(log10T9Root.Value), rho, Ye);
    }
  }

  return res;
}

NSEResult NSE::CalcFromEntropyAndDensityWithGuess(const double entropy,
    const double rho, const double Ye, const double T9Guess,
    const MuLM muLMGuess) {
  return CalcNSE(NSEMode::FixedRho, T9Guess, rho, Ye, entropy, muLMGuess);
}

NSEResult NSE::CalcFromInternalEnergyAndDensity(const double internalEnergy,
    const double rho, const double Ye) {
  return CalcFromInternalEnergyAndDensityWithTemperatureGuess(internalEnergy,
      rho, Ye, 10.0);
}

NSEResult NSE::CalcFromInternalEnergyAndDensityWithTemperatureGuess(
    const double internalEnergy, const double rho, const double Ye,
    const double T9Guess) {
  // we only have the internal energy and density, so we have to do a bisection
  // to find the temperature that gives this internal energy for the NSE
  // abundances obtained with that temperature
  auto f = [&] (const double log10T9) {
    double T9 = exp10(log10T9);
    auto res = CalcFromTemperatureAndDensity(T9, rho, Ye);
    auto eosRes = mpEOS->FromTAndRho(T9, rho, res.Y(), mNuclib);
    return eosRes.U() - internalEnergy;
  };

  BisectionReturnValue log10T9Root(-1.0, false);
  try {
    log10T9Root = BisectionMethod::FindRoot(f, log10(T9Guess), log10(1.001),
        log10(199.99));
  } catch (std::exception& ex) {
    throw std::runtime_error("Bisection failed in "
        "NSE::CalcFromInternalEnergyAndDensityWithTemperatureGuess with "
        "internal energy = " + std::to_string(internalEnergy) + ", rho = "
        + std::to_string(rho) + ", Ye = " + std::to_string(Ye) + "\nCause: "
        + ex.what());
  }

  auto res = CalcFromTemperatureAndDensity(exp10(log10T9Root.Value), rho, Ye);

  return NSEResult(res.Y(), exp10(log10T9Root.Value), rho, res.NumIterations(),
      res.Error(), res.GetMuLM());
}

MuLM NSE::CalcInitialGuess(const double T9, const double rho,
      const double Ye) {
  CheckInput(T9, rho, Ye);

  int Zl = mOpts.BasisZN0[0];
  int Nl = mOpts.BasisZN0[1];
  int Zm = mOpts.BasisZN1[0];
  int Nm = mOpts.BasisZN1[1];
  double fac = (double)(Zm * Nl - Zl * Nm);

  MuLM guess;
  guess.muL = 0.0;

  double TMeV = Constants::BoltzmannConstantInMeVPerGK * T9;
  double tempDensityFactor = GetTempDensityFactor(TMeV, rho);

  // find most bound nuclide
  double mostBindingEnergyPerNucleon = -0.1;
  int idxMostBound = -1;
  for (int i = 0; i < mNuclib.NumNuclides(); ++i) {
    //printf("%i: mass = %.10E\n", i, nuclideLibrary.MassesInMeV()[i]);
    double bindingEnergyPerNucleon = mNuclib.BindingEnergiesInMeV()[i]
        / (double)mNuclib.As()[i];
    if (bindingEnergyPerNucleon > mostBindingEnergyPerNucleon) {
      mostBindingEnergyPerNucleon = bindingEnergyPerNucleon;
      idxMostBound = i;
    }
  }

  //printf("Most bound nuclide: %s\n",
  //    nuclideLibrary.Names()[idxMostBound].c_str());

  // suppose the mass fraction of the most bound nuclide was 1
  double muMostBound = -mNuclib.BindingEnergiesInMeV()[idxMostBound]
      - TMeV * log((double)mNuclib.As()[idxMostBound]
          * mNuclib.GroundStateSpinTerms()[idxMostBound]
          * mNuclib.PartitionFunctionsWithoutSpinTerm(T9)[idxMostBound]
      * pow(mNuclib.MassesInMeV()[idxMostBound], 3.0 / 2.0)
      * tempDensityFactor);

  //printf("muMostBound = %.10E\n", muMostBound);

  int ZMostBound = mNuclib.Zs()[idxMostBound];
  int NMostBound = mNuclib.Ns()[idxMostBound];
  double muMFromMostBound = (muMostBound * fac
      - (double)(NMostBound * Zm - ZMostBound * Nm) * guess.muL)
      / (double)(ZMostBound * Nl - NMostBound * Zl);

  //printf("muMFromMostBound = %.10E\n", muMFromMostBound);

  // now suppose that the mass fractions of nucleons is 1, hence Y_p + Y_n = 1
  double Yp = Ye;
  double Yn = 1 - Ye;

  double muP = -TMeV
      * log(2.0 * pow(mNuclib.ProtonMassInMeV(), 3.0 / 2.0)
      * tempDensityFactor / Yp);
  double muN = -TMeV
      * log(2.0 * pow(mNuclib.NeutronMassInMeV(), 3.0 / 2.0)
      * tempDensityFactor / Yn);

  //printf("tempDensityFactor = %.10E, muP = %.10E, muN = %.10E\n",
  //    tempDensityFactor, muP, muN);

  double muMFromNucleons = (double)Zm * muP + (double)Nm * muN;

  // guess for muM is the larger of the two possibilities
  guess.muM = std::min(muMFromMostBound, muMFromNucleons);

  return guess;
}

void NSE::CalcY(const double T9, const double etaL, const double etaM,
    const double rho, const int n, const std::vector<double>& alpha,
    const std::vector<double>& beta, std::vector<double> * const pY,
    const std::vector<double> * const pScreenCorr) {
  double TMeV = Constants::BoltzmannConstantInMeVPerGK * T9;
  double tempDensityFactor = GetTempDensityFactor(TMeV, rho);
  auto gs = mNuclib.PartitionFunctionsWithoutSpinTerm(T9);

  if ((pScreenCorr == nullptr) || (!mOpts.DoScreening)) {
    for (int i = 0; i < n; ++i) {
      //printf("%i ", i);
      double eta = alpha[i] * etaL + beta[i] * etaM;
      //printf("TMeV = %.20E, eta = %.20E\n", TMeV, eta);
      double expArg = eta + mNuclib.BindingEnergiesInMeV()[i] / TMeV;
      //expArgs[i] = std::min(expArg, 100.0);
      expArg = std::min(expArg, 100.0);

      double m = mNuclib.MassesInMeV()[i];
      double G = mNuclib.GroundStateSpinTerms()[i] * gs[i];
      (*pY)[i] = exp(expArg) * G * pow(m, 3.0 / 2.0) * tempDensityFactor;
    }
  } else {
    for (int i = 0; i < n; ++i) {
      double eta = alpha[i] * etaL + beta[i] * etaM;
      double expArg = eta + mNuclib.BindingEnergiesInMeV()[i] / TMeV
          - (*pScreenCorr)[mNuclib.Zs()[i]];
      expArg = std::min(expArg, 100.0);

      double m = mNuclib.MassesInMeV()[i];
      double G = mNuclib.GroundStateSpinTerms()[i] * gs[i];
      (*pY)[i] = exp(expArg) * G * pow(m, 3.0 / 2.0) * tempDensityFactor;
    }
  }

  for (int i = 0; i < n; ++i) {
    if ((*pY)[i] > 1.0)
      (*pY)[i] = 1.0;
  }
}

double NSE::CalcF(const NSEMode mode, const double T9, const double rho,
    const double Ye, const double entropy, const int n,
    const std::vector<double>& Y, std::vector<double> * const pFs) {
  (*pFs)[0] = 1.0; // fA
  (*pFs)[1] = Ye; // fZ

  for (int i = 0; i < n; ++i) {
    (*pFs)[0] -= (double)mNuclib.As()[i] * Y[i];
    (*pFs)[1] -= (double)mNuclib.Zs()[i] * Y[i];
  }

  double f = (*pFs)[0] * (*pFs)[0] + (*pFs)[1] * (*pFs)[1];

  if (mode != NSEMode::FixedT9AndRho) {
    auto res = mpEOS->FromTAndRho(T9, rho, Y, mNuclib);
    (*pFs)[2] = res.S() / entropy - 1.0; // fS

    f += (*pFs)[2] * (*pFs)[2];
    return f / 3.0;
  } else {
    return 0.5 * f;
  }
}

void NSE::CalcJacobian(const NSEMode mode, const double T9,
    const double etaL, const double etaM, const double rho,
    const double entropy, const int n, const std::vector<double>& alpha,
    const std::vector<double>& beta, const std::vector<double>& Y,
    std::vector<std::vector<double>> * const pJ) {
  (*pJ)[0] = std::vector<double>((*pJ)[0].size(), 0.0);
  (*pJ)[1] = std::vector<double>((*pJ)[1].size(), 0.0);

  // compute analytic partial derivatives (not including screening)
  for (int i = 0; i < n; ++i) {
    double AY = (double)mNuclib.As()[i] * Y[i];
    double ZY = (double)mNuclib.Zs()[i] * Y[i];

    (*pJ)[0][0] -= AY * alpha[i];
    (*pJ)[0][1] -= AY * beta[i];
    (*pJ)[1][0] -= ZY * alpha[i];
    (*pJ)[1][1] -= ZY * beta[i];
  }

  if (mode == NSEMode::FixedRho) {
    double TMeV = Constants::BoltzmannConstantInMeVPerGK * T9;
    auto dLogGByLogTs =
        mNuclib.LogDerivativeOfPartitionFunctionsWithSpinTerm(T9); // TODO add option for log10Cap

    std::vector<double> dYiDT(n);
    std::vector<double> dsDYi(n);

    double T9Inverse = 1.0 / T9;
    for (int i = 0; i < n; ++i) {
      double eta = alpha[i] * etaL + beta[i] * etaM;
      double expArg = eta + mNuclib.BindingEnergiesInMeV()[i] / TMeV;
      expArg = std::min(expArg, 100.0);
      dYiDT[i] = Y[i] * T9Inverse
          * (1.5 - (expArg - eta) + dLogGByLogTs[i]);
      dsDYi[i] = 1.5 - expArg + dLogGByLogTs[i];
    }

    (*pJ)[2] = std::vector<double>((*pJ)[1].size(), 0.0);

    auto res = mpEOS->FromTAndRho(T9, rho, Y, mNuclib);

    (*pJ)[2][2] = res.dSdT9();

    for (int i = 0; i < n; ++i) {
      (*pJ)[0][2] -= (double)mNuclib.As()[i] * dYiDT[i];
      (*pJ)[1][2] -= (double)mNuclib.Zs()[i] * dYiDT[i];
      (*pJ)[2][0] += dsDYi[i] * Y[i] * alpha[i];
      (*pJ)[2][1] += dsDYi[i] * Y[i] * beta[i];
      (*pJ)[2][2] += dsDYi[i] * dYiDT[i];
    }

    (*pJ)[2][0] /= entropy;
    (*pJ)[2][1] /= entropy;
    (*pJ)[2][2] /= entropy;
  }
}

NSEResult NSE::DoNRIteration(const NSEMode mode, const double inputT9,
    const double inputRho, const double Ye, const double entropy, const int n,
    const std::vector<double>& alpha, const std::vector<double>& beta,
    MuLM initialMuLM, const int maxIterations,
    const std::vector<double> * const pScreenCorr) {
  std::vector<double> fs;
  std::vector<std::vector<double>> j;

  if (mode == NSEMode::FixedT9AndRho) {
    fs = std::vector<double>(2);
    j.push_back(std::vector<double>(2));
    j.push_back(std::vector<double>(2));
  } else {
    fs = std::vector<double>(3);
    j.push_back(std::vector<double>(3));
    j.push_back(std::vector<double>(3));
    j.push_back(std::vector<double>(3));
  }

  double T9OrRho;
  if (mode == NSEMode::FixedT9AndRho)
    T9OrRho = 0.0;
  else if (mode == NSEMode::FixedT9)
    T9OrRho = inputRho;
  else if (mode == NSEMode::FixedRho)
    T9OrRho = inputT9;
  else
    throw std::runtime_error("Invalid NSEMode");

  int iter = 0;
  std::vector<double> Y(n);

  //initial guess
  double etaL = initialMuLM.muL
      / (Constants::BoltzmannConstantInMeVPerGK * inputT9);
  double etaM = initialMuLM.muM
      / (Constants::BoltzmannConstantInMeVPerGK * inputT9);

  while (iter < maxIterations) {
    auto t9Rho = GetT9Rho(mode, inputT9, inputRho, T9OrRho);

    CalcY(t9Rho.T9, etaL, etaM, t9Rho.Rho, n, alpha, beta, &Y, pScreenCorr);
    CalcF(mode, t9Rho.T9, t9Rho.Rho, Ye, entropy, n, Y, &fs);

//    printf("etaL = %.5e, etaM = %.5e, fA = %.5e, fZ = %.5e\n", etaL, etaM,
//        fs[0], fs[1]);

    // check if we're converged
    double error = std::max(fabs(fs[0]), fabs(fs[1]));
    if (mode != NSEMode::FixedT9AndRho)
      error = std::max(error, fabs(fs[2]));

    if (error < mOpts.ConvergenceThreshold) {
      MuLM muLM;
      muLM.muL = etaL * Constants::BoltzmannConstantInMeVPerGK * t9Rho.T9;
      muLM.muM = etaM * Constants::BoltzmannConstantInMeVPerGK * t9Rho.T9;

      return NSEResult(Y, t9Rho.T9, t9Rho.Rho, iter, error, muLM);
    }

    CalcJacobian(mode, t9Rho.T9, etaL, etaM, t9Rho.Rho, entropy, n, alpha,
        beta, Y, &j);

//    printf(
//        "%10.4E  %10.4E\n%10.4E  %10.4E\n",
//        j[0][0], j[0][1],
//        j[1][0], j[1][1]);

    double deltaEtaL = 0.0;
    double deltaEtaM = 0.0;
    double deltaT9OrRho = 0.0;

    double det = 0.0;
    if (mode == NSEMode::FixedT9AndRho) {
      det = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    } else {
      det = Det(
          j[0][0], j[0][1], j[0][2],
          j[1][0], j[1][1], j[1][2],
          j[2][0], j[2][1], j[2][2]);
    }

//    printf("det J = %.5e\n\n", det);
    if (fabs(det) < 1.0E-100)
      throw std::runtime_error("Singular Jacobian in NSE");

    if (mode == NSEMode::FixedT9AndRho) {
//      printf(
//          "%10.4E  %10.4E\n%10.4E  %10.4E\n\n",
//          j[0][0], j[0][1],
//          j[1][0], j[1][1]);

      deltaEtaL = (j[1][1] * fs[0] - j[0][1] * fs[1]) / det;
      deltaEtaM = (-j[1][0] * fs[0] + j[0][0] * fs[1]) / det;
    } else if (mode == NSEMode::FixedRho) {
      deltaEtaL = Det(
          fs[0], j[0][1], j[0][2],
          fs[1], j[1][1], j[1][2],
          fs[2], j[2][1], j[2][2]) / det;

      deltaEtaM = Det(
          j[0][0], fs[0], j[0][2],
          j[1][0], fs[1], j[1][2],
          j[2][0], fs[2], j[2][2]) / det;

      deltaT9OrRho = Det(
          j[0][0], j[0][1], fs[0],
          j[1][0], j[1][1], fs[1],
          j[2][0], j[2][1], fs[2]) / det;

//      printf(
//          "%10.4E  %10.4E  %10.4E\n%10.4E  %10.4E  %10.4E\n%10.4E  %10.4E  %10.4E\n\n",
//          j[0][0], j[0][1], j[0][2],
//          j[1][0], j[1][1], j[1][2],
//          j[2][0], j[2][1], j[2][2]);
    } else {
      throw std::runtime_error("Invalid NSEMode");
    }

    // scale the update so that the change is not too large
    double fac = std::min(6.9 / (fabs(deltaEtaL) + 1.0E-20),
        6.9 / (fabs(deltaEtaM) + 1.0E-20));

    if (mode == NSEMode::FixedRho) {
      //printf("T9 = %.20E, deltaT9 = %.4E\n", T9OrRho, deltaT9OrRho);
      fac = std::min(fac, 0.1 * T9OrRho / (fabs(deltaT9OrRho) + 1.0E-20));
    }

    fac = std::min(std::max(fac, mOpts.StepSizeLimit), 1.0);

//    if (mode == NSEMode::FixedRho) {
//      printf("values = (%12.5E, %12.5E, %12.5E), error = (%12.5E, %12.5E, %12.5E), delta = (%12.5E, %12.5E, %12.5E), fac = %12.5E\n",
//          etaL, etaM, T9OrRho, fs[0], fs[1], fs[2], deltaEtaL, deltaEtaM, deltaT9OrRho, fac);
//    }

    etaL -= fac * deltaEtaL;
    etaM -= fac * deltaEtaM;
    T9OrRho -= fac * deltaT9OrRho;
    if (mode == NSEMode::FixedRho)
      T9OrRho = std::max(T9OrRho, 1.0); // don't let T9 go below 1 GK

    ++iter;
  }

  throw NSENotConverged();
}

NSEResult NSE::CalcNSE(const NSEMode mode, const double inputT9,
    const double inputRho, const double Ye, const double entropy,
    const MuLM muLMGuess) {
  if (mode == NSEMode::FixedT9AndRho) {
    if (!std::isnan(entropy))
      throw std::invalid_argument("In fixed T9 and rho mode, entropy must be "
          "NaN");
  } else {
    if (std::isnan(entropy))
      throw std::invalid_argument("Entropy can be NaN "
          "only in fixed T9 and rho mode");
  }

  if (mode == NSEMode::FixedT9)
    throw std::runtime_error("fixed T9 mode is not implemented yet");

  CheckInput(inputT9, inputRho, Ye);

  MuLM guess = muLMGuess;
  bool guessWasProvided = true;
  if (std::isnan(muLMGuess.muL) || std::isnan(muLMGuess.muM)) {
    guess = CalcInitialGuess(inputT9, inputRho, Ye);
    guessWasProvided = false;
  }

  // calculate alphas and betas
  int Zl = mOpts.BasisZN0[0];
  int Nl = mOpts.BasisZN0[1];
  int Zm = mOpts.BasisZN1[0];
  int Nm = mOpts.BasisZN1[1];
  double fac = 1.0 / (double)(Zl * Nm - Nl * Zm);

  int n = mNuclib.NumNuclides();
  std::vector<double> alpha(n);
  std::vector<double> beta(n);

  for (int i = 0; i < n; ++i) {
    int Zi = mNuclib.Zs()[i];
    int Ni = mNuclib.Ns()[i];
    alpha[i] = (double)(Zi * Nm - Ni * Zm) * fac;
    beta[i] = (double)(Zi * Nl - Ni * Zl) * (-fac);
  }

  int maxIterations =
      guessWasProvided ? mOpts.MaxIterationsWithGuess : mOpts.MaxIterations;

  try {
    if (mOpts.DoScreening) {
      // first get NSE without screening
      auto res = DoNRIteration(mode, inputT9, inputRho, Ye, entropy, n,
          alpha, beta, guess, maxIterations, nullptr);

      int totalIter = res.NumIterations();

      // now iterate computing screening corrections and calculating NSE with
      // the fixed screening corrections
      std::vector<double> prevY(res.Y());
      double error = 0.0;
      int iter = 0;

      do {
        auto therm = mpEOS->FromTAndRho(inputT9, inputRho, prevY, mNuclib);
        auto scr = mpScreening->Compute(therm, prevY, mNuclib);

        auto newRes = DoNRIteration(mode, inputT9, inputRho, Ye, entropy, n,
            alpha, beta, res.GetMuLM(), maxIterations, &scr.Mus());
        totalIter += newRes.NumIterations();

        error = 0.0;
        for (int i = 0; i < n; ++i) {
          error = std::max(error, fabs(prevY[i] - newRes.Y()[i]));
        }
//        printf("%i: screening error = %.14e\n", iter, error);

        if (error < mOpts.ScreeningConvergenceThreshold) {
          newRes.mNumIterations = totalIter;
          return newRes;
        }

        ++iter;
        prevY = newRes.Y();
        res = newRes;
      } while ((error >= mOpts.ScreeningConvergenceThreshold)
          && (iter < mOpts.MaxScreeningIterations));

      // failed, continue below
    } else {
      return DoNRIteration(mode, inputT9, inputRho, Ye, entropy, n, alpha, beta,
          guess, maxIterations, nullptr);
    }
  } catch (NSENotConverged&) {
    // failed, continue below
  }

  // if a guess was provided, try again without the guess
  if (guessWasProvided)
    return CalcNSE(mode, inputT9, inputRho, Ye, entropy, MuLM());

  std::string modeStr;
  if (mode == NSEMode::FixedT9AndRho) {
    modeStr = "FixedT9AndRho";
  } else if (mode == NSEMode::FixedT9) {
    modeStr = "FixedT9 (entropy = " + std::to_string(entropy) + ")";
  } else if (mode == NSEMode::FixedRho) {
    modeStr = "FixedRho (entropy = " + std::to_string(entropy) + ")";
  } else {
    throw std::runtime_error("Invalid NSEMode");
  }

  throw std::runtime_error("Reached maximum number of iterations in NSE. "
      "inputT9 = " + std::to_string(inputT9) + ", inputRho = "
      + std::to_string(inputRho) + ", Ye = " + std::to_string(Ye)
      + ", mode = " + modeStr);
}

void NSE::CheckInput(const double T9, const double rho, const double Ye) {
  if ((fcmp(T9, 200.0) > 0) || (fcmp(T9, 1.0) < 0))
    throw std::invalid_argument("T9 for NSE must be between 1.0 and 100.0, "
        "but got T9 = " + std::to_string(T9));
  if ((fcmp(rho, 1.0E13) > 0) || (fcmp(rho, 1.0) < 0))
    throw std::invalid_argument("rho for NSE must be between 1.0 and 1.0E13, "
        "but got rho = " + std::to_string(rho));
  if ((fcmp(Ye, 0.999) > 0) || (fcmp(Ye, 0.001) < 0))
    throw std::invalid_argument("Ye for NSE must be between 0.001 and 0.999, "
        "but got Ye = " + std::to_string(Ye));
}

