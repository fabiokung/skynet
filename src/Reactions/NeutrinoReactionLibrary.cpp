/// \file NeutrinoReactionLibrary.cpp
/// \author lroberts
/// \since May 12, 2015
///
/// \brief
///
///

#include "Reactions/NeutrinoReactionLibrary.hpp"
#include "Utilities/FunctionIntegrator.hpp"
#include "Utilities/Constants.hpp"
#include "Reactions/ReactionData.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>

#define MAX_EXPONENTIAL 200.0

namespace { // unnamed so this can only be used in this file

const double BaryonMassInMeV = 939.0;

double WeakMagnetismCorrection(const double eNuMeV, const bool nueCap) {
  // Charged-current recoil and weak-magnetism correction from
  // Horowitz, Phys. Rev. D 65, 043001 (2002), Eqs. (24) and (25).
  // These first-order factors are accurate at the supernova neutrino
  // energies used by the νp-process example.
  const double coeff = nueCap ? 1.01 : -7.22;
  return std::max(0.0, 1.0 + coeff * eNuMeV / BaryonMassInMeV);
}

double WeakMagnetismExactCorrection(const double eNuMeV, const bool nueCap) {
  // Charged-current recoil + weak-magnetism correction to all orders in
  // e = E_nu/M, Horowitz, Phys. Rev. D 65, 043001 (2002), Eq. (22). The
  // (1+2e)^3 denominator keeps R(k) positive at all energies, unlike the
  // first-order form which goes negative for ν̄_e above ~50 MeV.
  // Table I charged-current couplings: cV = 1, cA = g_A, F2 = κ_p - κ_n.
  const double cV = 1.0, cA = 1.26, F2 = 3.706;
  const double e = eNuMeV / BaryonMassInMeV;
  const double interf = nueCap ? 1.0 : -1.0; // ± axial-vector interference
  const double num = cV * cV * (1.0 + 4.0 * e + 16.0 / 3.0 * e * e)
      + 3.0 * cA * cA * (1.0 + 4.0 / 3.0 * e) * (1.0 + 4.0 / 3.0 * e)
      + interf * 4.0 * (cV + F2) * cA * e * (1.0 + 4.0 / 3.0 * e)
      + 8.0 / 3.0 * cV * F2 * e * e
      + 5.0 / 3.0 * e * e * (1.0 + 2.0 / 5.0 * e) * F2 * F2;
  const double den = (cV * cV + 3.0 * cA * cA) * pow(1.0 + 2.0 * e, 3);
  return num / den;
}

double CaptureCorrection(const double eNuMeV, const bool nueCap,
    const NeutrinoCorrectionMode correctionMode) {
  switch (correctionMode) {
    case NeutrinoCorrectionMode::WeakMagnetism:
      return WeakMagnetismCorrection(eNuMeV, nueCap);
    case NeutrinoCorrectionMode::WeakMagnetismExact:
      return WeakMagnetismExactCorrection(eNuMeV, nueCap);
    case NeutrinoCorrectionMode::None:
    default:
      return 1.0;
  }
}

std::function<double(double)> MakeECapFunc(const double scaledT,
    const double eta, const double meScaled, const double scaledQ,
    const int nuPow, const double scaledWm, const double eScale,
    const bool nueCap, const NeutrinoCorrectionMode correctionMode,
    const std::function<double(double)> distributionFunction) {
    return [=] (const double x) {
        double ae = x / scaledT - eta;
        if (ae > MAX_EXPONENTIAL)
          return 0.0;
        ae = exp(ae);
        const double eNuMeV = (x - scaledQ) * eScale;
        const double correction = CaptureCorrection(eNuMeV, nueCap,
            correctionMode);

        return x * sqrt(x * x - meScaled * meScaled)
            * pow((x - scaledQ), nuPow)
            * (1.0 + x * scaledWm) * correction / (1.0 + ae)
            * (1.0 - distributionFunction(eNuMeV));
      };
}

std::function<double(double)> MakeNuCapFunc(const double scaledT,
    const double eta, const double meScaled, const double scaledQ,
    const int nuPow, const double scaledWm, const double eScale,
    const bool nueCap, const NeutrinoCorrectionMode correctionMode,
    const std::function<double(double)> distributionFunction) {
    return [=] (const double x) {
        double ae = x / scaledT - eta;
        if (ae > MAX_EXPONENTIAL)
          ae = MAX_EXPONENTIAL;
        ae = exp(ae);
        const double eNuMeV = (x - scaledQ) * eScale;
        const double correction = CaptureCorrection(eNuMeV, nueCap,
            correctionMode);

        return x * sqrt(x * x - meScaled * meScaled)
            * pow((x - scaledQ), nuPow)
            * (1.0 + x * scaledWm) * correction * ae / (1.0 + ae)
            * distributionFunction(eNuMeV);
      };
}

std::function<double(double)> MakeDecayFunc(const double scaledT,
    const double eta, const double meScaled, const double scaledQ,
    const int nuPow, const double scaledWm, const double eScale,
    const bool nueCap, const NeutrinoCorrectionMode correctionMode,
    const std::function<double(double)> distributionFunction) {
    return [=] (const double x) {
        double ae = x / scaledT - eta;
        if (ae > MAX_EXPONENTIAL)
          ae = MAX_EXPONENTIAL;
        ae = exp(ae);
        const double eNuMeV = (x - scaledQ) * eScale;
        const double correction = CaptureCorrection(eNuMeV, nueCap,
            correctionMode);

        return x * sqrt(x * x - meScaled * meScaled)
            * pow((x + scaledQ), nuPow)
            * (1.0 + x * scaledWm) * correction * ae / (1.0 + ae)
            * (1.0 - distributionFunction(eNuMeV));
      };
}
} // namespace [unnamed]

NeutrinoReactionLibrary::NeutrinoReactionLibrary(const Neutrino neutrinoLib,
    const std::string& description,
    const NuclideLibrary& nucLib, const NetworkOptions& opts, bool onlyNuCap,
    bool nuHeating, bool includeBeta, NeutrinoCorrectionMode correctionMode) :
    ReactionLibraryBase(ReactionType::Weak, description,
        neutrinoLib.GetSource(), neutrinoLib.GetValidReactions(nucLib),
        nucLib, opts, false),
    mOnlyNuCap(onlyNuCap),
    mNuHeating(nuHeating),
    mIncludeBeta(includeBeta),
    mCorrectionMode(correctionMode) {

  auto entries = neutrinoLib.Entries();
  std::vector<double> Q, matrixElement, Wm;
  std::vector<bool> nueCap;
  for (auto entry : entries) {
    Q.push_back(entry.GetQ());
    matrixElement.push_back(entry.GetMatrixElement());
    Wm.push_back(entry.GetWm());
    nueCap.push_back(entry.IsNueReaction());
  }

  if (mCorrectionMode != NeutrinoCorrectionMode::None) {
    for (auto wm : Wm) {
      if (wm != 0.0) {
        throw std::runtime_error("Analytic weak-magnetism correction modes "
            "cannot be combined with non-zero neutrino reaction "
            "weak-magnetism coefficients.");
      }
    }
  }

  mQ = ReactionData<double>(Q);
  mMatrixElement = ReactionData<double>(matrixElement);
  mWm = ReactionData<double>(Wm);
  mNueCap = ReactionData<bool>(nueCap);

  mRates = std::vector<double>(NumAllReactions());
  mInverseRates = std::vector<double>(NumAllReactions());
  mHeatingRates = std::vector<double>(NumAllReactions());
  mHeatingInverseRates = std::vector<double>(NumAllReactions());
}

void NeutrinoReactionLibrary::PrintAdditionalInfo(
    NetworkOutput * const pOutput) const {
  pOutput->Log("#   Only nu capture: %s\n", mOnlyNuCap ? "yes" : "no");
  pOutput->Log("#   Compute heating: %s\n", mNuHeating ? "yes" : "no");
  const char* modeName = "None";
  if (mCorrectionMode == NeutrinoCorrectionMode::WeakMagnetism)
    modeName = "WeakMagnetism";
  else if (mCorrectionMode == NeutrinoCorrectionMode::WeakMagnetismExact)
    modeName = "WeakMagnetismExact";
  pOutput->Log("#   Correction mode: %s\n", modeName);
}

void NeutrinoReactionLibrary::DoLoopOverReactionData(
    const std::function<void(GeneralReactionData * const)>& func) {
  func(&mQ);
  func(&mMatrixElement);
  func(&mWm);
  func(&mNueCap);
  mRates = std::vector<double>(Reactions().ActiveData().size());
  mInverseRates = std::vector<double>(Reactions().ActiveData().size());
  mHeatingRates = std::vector<double>(Reactions().ActiveData().size());
  mHeatingInverseRates = std::vector<double>(Reactions().ActiveData().size());
}

void NeutrinoReactionLibrary::DoCalculateRates(
    const ThermodynamicState thermoState,
    const std::vector<double>& /*partitionFunctionsWithoutSpinTerms*/,
    const double /*expArgumentCap*/, const std::vector<int> * const /*pZs*/,
    const std::vector<double> * const /*pScreeningChemicalPotentialCorrection*/) {

  if (thermoState.T9() < Options().MinT9ForNeutrinoReactions) {
    std::fill(mRates.begin(), mRates.end(), 0.0);
    std::fill(mInverseRates.begin(), mInverseRates.end(), 0.0);
    return;
  }

  FunctionIntegrator integrator;

  // Define the forward function
  //double eScale    = std::max(Constants::ElectronMassInMeV,
  //    thermoState.T9()*Constants::BoltzmannConstantInMeVPerGK);
  double eScale = Constants::ElectronMassInMeV;
  double scaledQ = 0.0;
  double scaledWm = 0.0;
  double scaledT = thermoState.T9() * Constants::BoltzmannConstantInMeVPerGK
      / eScale;
  double meScaled = Constants::ElectronMassInMeV / eScale;
  // Electron chemical degeneracy parameter corrected for the electron rest mass
  double etae = thermoState.EtaElectron() + Constants::ElectronMassInMeV
      / (thermoState.T9() * Constants::BoltzmannConstantInMeVPerGK);

  std::shared_ptr<NeutrinoDistribution> nuDist =
      thermoState.GetNeutrinoDistribution();

  // Calculate the neutrino "luminosity" for consistency check
  // Should deviate from the actual luminosity since the quantity
  // we are calculating is proportional to the energy density, but
  // the deviation should be very small when the radius is much
  // greater than the neutrino sphere radius
  //auto lum_func = [&nuDist, &nuSpec] (const double& x) {
  //  double f =  nuDist->GetDistributionFunction(x,nuSpec);
  //  return x*x*x*f;
  //};
  //double Lnu = integrator.Integrate(lum_func,0.0, 300.0);
  //Lnu = Lnu/pow(197.3,3)*1.e39*1.602e-6*3.e10
  //  * 2.0/3.14159*nuDist->GetRadius()*nuDist->GetRadius();
  //std::cout << "Lnu " << Lnu  << " " << nuDist->GetRadius() << " "
  //    << nuDist->GetRnu(nuSpec) << std::endl;

  const double rate_const = pow(eScale / Constants::ElectronMassInMeV, 5)
      * log(2.0) / mK;
  const double heat_const = rate_const * eScale; // MeV/s

  for (unsigned int i = 0; i < mQ.size(); i++) {
    scaledQ = mQ[i] / eScale;
    scaledWm = mWm[i] * eScale;
    double lower_lim = std::max(Constants::ElectronMassInMeV / eScale, scaledQ);

    double eta = etae;
    NeutrinoSpecies nuSpec = NeutrinoSpecies::NuE;

    if (!mNueCap[i]) {
      nuSpec = NeutrinoSpecies::AntiNuE;
      eta = -etae;
    }

    double ecapInt = 0.0;
    double ecapHeatInt = 0.0;
    if (!mOnlyNuCap) {
      try {
        auto eCapFunc = MakeECapFunc(scaledT, eta, meScaled, scaledQ, 2,
            scaledWm, eScale, mNueCap[i], mCorrectionMode,
            nuDist->DistributionFunction(nuSpec));
        if (eta > -100.0) 
          ecapInt = integrator.Integrate(eCapFunc, lower_lim,
              std::numeric_limits<double>::infinity());
        if (-mQ[i] > lower_lim && mIncludeBeta) { 
          auto decayFunc = MakeDecayFunc(scaledT, eta, meScaled, scaledQ, 2,
              scaledWm, eScale, mNueCap[i], mCorrectionMode,
              nuDist->DistributionFunction(nuSpec));
          ecapInt = integrator.Integrate(decayFunc, lower_lim, -scaledQ);
        }
        if (mNuHeating) {
          auto eCapFunc = MakeECapFunc(scaledT, eta, meScaled, scaledQ, 3,
              scaledWm, eScale, mNueCap[i], mCorrectionMode,
              nuDist->DistributionFunction(nuSpec));
          if (eta > -100.0)
            ecapHeatInt = integrator.Integrate(eCapFunc, lower_lim,
                std::numeric_limits<double>::infinity());
          if (-mQ[i] > Constants::ElectronMassInMeV && mIncludeBeta) { 
            auto BetaFunc = MakeDecayFunc(scaledT, eta, meScaled, scaledQ, 3,
                scaledWm, eScale, mNueCap[i], mCorrectionMode,
                nuDist->DistributionFunction(nuSpec));
            ecapHeatInt += integrator.Integrate(BetaFunc, 
                Constants::ElectronMassInMeV/eScale,
                -scaledQ);
          }
        } else {
          ecapHeatInt = 0.0;
        }
      } catch (int e) {
        ecapInt = 0.0;
        ecapHeatInt = 0.0;
        std::cerr << "Electron capture integration error " << e <<  
            Reactions()[i].String() << std::endl;
      } catch (...) {
        ecapInt = 0.0;
        ecapHeatInt = 0.0;
        std::cerr << "Electron capture integration error " << 
            Reactions()[i].String() << std::endl;
      }
      mRates[i] = mMatrixElement[i] * rate_const * ecapInt;
      mHeatingRates[i] = mMatrixElement[i] * heat_const * ecapHeatInt;
    } else {
      mRates[i] = 0.0;
      mHeatingRates[i] = 0.0;
    }

    double nucapInt = 0.0;
    double nucapHeatInt = 0.0;
    try {
      auto nuCapFunc = MakeNuCapFunc(scaledT, eta, meScaled, scaledQ, 2,
          scaledWm, eScale, mNueCap[i], mCorrectionMode,
          nuDist->DistributionFunction(nuSpec));
      nucapInt = integrator.Integrate(nuCapFunc, lower_lim,
          std::numeric_limits<double>::infinity());
      if (mNuHeating) {
        auto nuCapFunc = MakeNuCapFunc(scaledT, eta, meScaled, scaledQ, 3,
            scaledWm, eScale, mNueCap[i], mCorrectionMode,
            nuDist->DistributionFunction(nuSpec));
        nucapHeatInt = -integrator.Integrate(nuCapFunc, lower_lim,
            std::numeric_limits<double>::infinity());
      } else {
        nucapHeatInt = 0.0;
      }
    } catch (int e) {
      nucapInt = 0.0;
      nucapHeatInt = 0.0;
      std::cerr << "Neutrino capture integration error " << e << std::endl;
    } catch (...) {
      nucapInt = 0.0;
      nucapHeatInt = 0.0;
      std::cerr << "Neutrino capture integration error " << std::endl;
    }
    mInverseRates[i] = mMatrixElement[i] * rate_const * nucapInt;
    mHeatingInverseRates[i] = mMatrixElement[i] * heat_const * nucapHeatInt;
  }
}
