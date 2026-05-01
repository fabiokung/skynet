/// \file ReadFFN.cpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#include <vector>

#include "BuildInfo.hpp"
#include "Utilities/FloatingPointExceptions.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Neutrino.hpp"
#include "Reactions/NeutrinoReactionLibrary.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NSEOptions.hpp"
#include "Utilities/FunctionVsTimeWrapper.hpp"
#include "Utilities/FloatingPointExceptions.hpp"
#include "Utilities/Constants.hpp"

#include <iostream>
#include <cmath>

namespace {

class NarrowNeutrinoDistribution : public NeutrinoDistribution {
public:
  NarrowNeutrinoDistribution(const double centerMeV, const double widthMeV,
      const double normalization) :
      NeutrinoDistribution({NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE}),
      mCenterMeV(centerMeV),
      mWidthMeV(widthMeV),
      mNormalization(normalization) {}

  std::function<double(double)> DistributionFunction(
      const NeutrinoSpecies /*species*/) const {
    return [=] (const double eNuMeV) {
      const double x = (eNuMeV - mCenterMeV) / mWidthMeV;
      return mNormalization * std::exp(-x * x);
    };
  }

private:
  double mCenterMeV;
  double mWidthMeV;
  double mNormalization;
};

class NarrowNeutrinoHistory : public NeutrinoHistory {
public:
  NarrowNeutrinoHistory(const double centerMeV, const double widthMeV,
      const double normalization) :
      NeutrinoHistory(true, std::vector<NeutrinoSpecies>{
          NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE}),
      mCenterMeV(centerMeV),
      mWidthMeV(widthMeV),
      mNormalization(normalization) {}

  std::shared_ptr<NeutrinoDistribution> operator()(
      const double /*time*/) const {
    return std::shared_ptr<NeutrinoDistribution>(
        new NarrowNeutrinoDistribution(mCenterMeV, mWidthMeV,
            mNormalization));
  }

  std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const {
    return
        std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>(
            new NarrowNeutrinoHistory(*this));
  }

  std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const {
    return std::shared_ptr<NeutrinoHistory>(new NarrowNeutrinoHistory(*this));
  }

  void PrintInfo(NetworkOutput * const pOutput) const {
    pOutput->Log("# Neutrino History: narrow test spectrum\n");
  }

private:
  double mCenterMeV;
  double mWidthMeV;
  double mNormalization;
};

} // namespace

int main(int, char**) {

  FloatingPointExceptions::Enable();

  NetworkOptions opts;

  // Read in EoS table
  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");
  const double deltaNP = 1.293;

  // Set up neutron proton nuclear chart
  std::array<double, 24> partitionFunction;
  partitionFunction.fill(1.0);
  std::vector<Nuclide> nuclides = {
      Nuclide(0, 1, 0.0, 0.5, partitionFunction, "n"),
      Nuclide(1, 1, -deltaNP, 0.5, partitionFunction, "p")
  };
  NuclideLibrary nuclib(nuclides, "neutron proton");

  // Test building neutrino interactions from beta decay
  {
    // These are a couple of different ways of building electron captures
    // from beta decays

    //REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
    //    ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
    //    "Weak reactions", nuclib, opts);
    //
    //std::vector<NeutrinoEntry> nuRates;
    //for (int i=0; i<weakReactionLibrary.Reactions().size(); ++i) {
    //  std::cout << weakReactionLibrary.Reactions()[i].String() << " "
    //      << exp(weakReactionLibrary.RateFittingCoefficients()[0][i]) << std::endl;
    //  nuRates.push_back(NeutrinoEntry::FromVacuumBetaDecay(
    //      weakReactionLibrary.Reactions()[i],
    //      exp(weakReactionLibrary.RateFittingCoefficients()[0][i]),
    //      nuclib, false));
    //  nuRates.push_back(NeutrinoEntry::FromVacuumBetaDecay(
    //      weakReactionLibrary.Reactions()[i],
    //      exp(weakReactionLibrary.RateFittingCoefficients()[0][i]),
    //      nuclib, true));
    //}
    //NeutrinoReactionLibrary reactionLib(Neutrino(nuRates, "Test"),
    //    "Neutrino Reactions", nuclib, opts, false, false, true);

    //NeutrinoEntry NDecay = NeutrinoEntry::FromVacuumBetaDecay(
    //  Reaction({"n"},{"p"},{1},{1}, true, false, false, "Neutron Decay", nuclib),
    //  log(2.0),
    //  nuclib, false);
    //NeutrinoEntry PDecay = NeutrinoEntry::FromVacuumBetaDecay(
    //  Reaction({"n"},{"p"},{1},{1}, true, false, false, "Neutron Decay", nuclib),
    //  log(2.0),
    //  nuclib, true);
    //NeutrinoReactionLibrary reactionLib(Neutrino({NDecay, PDecay}, "Test"),
    //    "Neutrino Reactions", nuclib, opts, false, false, true);

    Neutrino nuRates = Neutrino::FromREACLIBDecays(SkyNetRoot + "/data/reaclib",
        nuclib, opts);
    NeutrinoReactionLibrary reactionLib(nuRates,
        "Neutrino Reactions", nuclib, opts, false, false, true);

    std::vector<double> yInit { 1.0-1.e-12, 1.e-12 };
    double TGK = 5.e-1;
    double RHOGCC = 1.e1;
    opts.MaxDt = 1.0;
    ReactionNetwork network(nuclib, { &reactionLib }, &helm, nullptr, opts);
    std::vector<double> yFinal =
        network.Evolve(yInit, 0.0, 1222.0, FVT([&TGK] (double) {return TGK;}),
            FVT([&RHOGCC] (double) {return RHOGCC;}), "SkyNet_output",
            1.e-10).FinalY();
    if (fabs(yFinal[0] - 0.25) > 0.01) return 1;
    opts = NetworkOptions();
  }

  // Read in neutrino reactions and make library
  Neutrino lib(SkyNetRoot + "/data/neutrino_reactions.dat", nuclib);
  std::vector<Reaction> reacs = lib.Reactions(nuclib);

  for (auto reac : reacs) {
    printf("%s\n", reac.Label().c_str());
    printf("%s\n", reac.String().c_str());
  }

  NeutrinoReactionLibrary reactionLib(lib, "Neutrino Reactions", nuclib, opts);

  // Setup constant background and neutrino conditions
  const double TGK = 5.e-1;
  const double RHOGCC = 1.e6;
  const double ENUE = 12.0 / 3.59714 / Constants::BoltzmannConstantInMeVPerGK;
  const double ENUB = 12.0 / 3.59714 / Constants::BoltzmannConstantInMeVPerGK;
  const double LNUE = 1.e52;
  const double LNUB = 1.e52;
  const double RAD = 1.e8;
  ReactionNetwork network(nuclib, { &reactionLib }, &helm, nullptr, opts);
  auto nuHist = NeutrinoHistoryBlackBody::CreateConstant({ 0.0, 1.e20 },
      { RAD, RAD }, { NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE },
      { ENUE, ENUB }, { 0.0, 0.0 }, { LNUE, LNUB }, false);

  {
    NeutrinoReactionLibrary defaultLib(lib, "Default Neutrino Reactions",
        nuclib, opts, true);
    NeutrinoReactionLibrary noneLib(lib, "Uncorrected Neutrino Reactions",
        nuclib, opts, true, false, false, NeutrinoCorrectionMode::None);
    NeutrinoReactionLibrary wmLib(lib, "Weak Magnetism Neutrino Reactions",
        nuclib, opts, true, false, false,
        NeutrinoCorrectionMode::WeakMagnetism);

    ThermodynamicState thermoState(TGK, RHOGCC, 0.0, 0.0, 0.0, 0.5, 0.0,
        nuHist(0.0));
    std::vector<double> partitionFunctions(nuclib.NumNuclides(), 1.0);

    defaultLib.CalculateRates(thermoState, partitionFunctions,
        opts.RateExpArgumentCap, nullptr, nullptr);
    noneLib.CalculateRates(thermoState, partitionFunctions,
        opts.RateExpArgumentCap, nullptr, nullptr);
    wmLib.CalculateRates(thermoState, partitionFunctions,
        opts.RateExpArgumentCap, nullptr, nullptr);

    for (unsigned int i = 0; i < defaultLib.InverseRates().size(); ++i) {
      if (fabs(defaultLib.InverseRates()[i] - noneLib.InverseRates()[i])
          > 1.e-12 * noneLib.InverseRates()[i])
        return EXIT_FAILURE;

      // Consistency check with Horowitz, Phys. Rev. D 65, 043001 (2002).
      // Weak magnetism enhances neutrino capture rates and suppresses
      // antineutrino capture rates.
      if (lib.Entries()[i].IsNueReaction()) {
        if (wmLib.InverseRates()[i] <= noneLib.InverseRates()[i])
          return EXIT_FAILURE;
      } else {
        if (wmLib.InverseRates()[i] >= noneLib.InverseRates()[i])
          return EXIT_FAILURE;
      }
    }

    auto narrow20MeVDist = std::shared_ptr<NeutrinoDistribution>(
        new NarrowNeutrinoDistribution(20.0, 2.0, 1.0));
    thermoState.SetNeutrinoDistribution(narrow20MeVDist);
    noneLib.CalculateRates(thermoState, partitionFunctions,
        opts.RateExpArgumentCap, nullptr, nullptr);
    wmLib.CalculateRates(thermoState, partitionFunctions,
        opts.RateExpArgumentCap, nullptr, nullptr);

    for (unsigned int i = 0; i < wmLib.InverseRates().size(); ++i) {
      if (lib.Entries()[i].IsNueReaction())
        continue;

      const double uncorrectedRate = noneLib.InverseRates()[i];
      if (uncorrectedRate <= 0.0) {
        std::cerr << "WeakMagnetism Horowitz qualitative check failed: "
            << "invalid uncorrected anti-neutrino capture rate for the 20 MeV "
            << "comparison." << std::endl;
        return EXIT_FAILURE;
      }
      const double rateRatio = wmLib.InverseRates()[i] / uncorrectedRate;
      if (!std::isfinite(rateRatio)) {
        std::cerr << "WeakMagnetism Horowitz qualitative check failed: "
            << "invalid anti-neutrino capture rate ratio for the 20 MeV "
            << "comparison." << std::endl;
        return EXIT_FAILURE;
      }
      const double reduction = 1.0 - rateRatio;
      // Qualitative check against Horowitz, Phys. Rev. D 65, 043001 (2002):
      // Eq. (24) gives R_nuebar ~= 1 - 7.22 k/m, and the text notes a
      // roughly 15% charged-current antineutrino opacity reduction at
      // k ~= 20 MeV. With m = 939 MeV this is 15.38%.
      const double horowitz20MeVReduction = 7.22 * 20.0 / 939.0;
      if (fabs(reduction - horowitz20MeVReduction) > 5.0e-3) {
        std::cerr << "WeakMagnetism Horowitz qualitative check failed: "
            << "expected about 15.38% anti-neutrino reduction at k=20 MeV, "
            << "got " << reduction * 100.0 << "%." << std::endl;
        return EXIT_FAILURE;
      }
    }
  }

  network.LoadNeutrinoHistory(nuHist.MakeSharedPtr());

  // Evolve the network
  std::vector<double> yInit { 0.1, 0.9 };
  std::vector<double> yFinal =
      network.Evolve(yInit, 0.0, 1.e10, FVT([&TGK] (double) {return TGK;}),
          FVT([&RHOGCC] (double) {return RHOGCC;}), "SkyNet_output",
          1.e-10).FinalY();

  // Analytically calculate what the equilibrium should be for neutrinos only
  const double trescale = 3.59714;
  const double nuePhaseSpace = LNUE * (23.3309 * pow(ENUE / trescale, 4)
      + 11.3644 * pow(ENUE / trescale, 3) * deltaNP
      + 1.80309 * pow(ENUE / trescale, 2)) / pow(ENUE, 3);
  const double nubPhaseSpace = LNUB * (23.3309 * pow(ENUB / trescale, 4)
      - 11.3644 * pow(ENUB / trescale, 3) * deltaNP
      + 1.80309 * pow(ENUB / trescale, 2)) / pow(ENUB, 3);
  const double approxYeFin = 1.0 / (1.0 + nubPhaseSpace / nuePhaseSpace);

  const double relativeErr = (approxYeFin - yFinal[1]) / yFinal[1];
  std::cout << yFinal[1] << " " << approxYeFin << " "
      << relativeErr << std::endl;

  if (relativeErr > 1.e-3) return EXIT_FAILURE;

  NeutrinoReactionLibrary narrowWmReactionLib(lib,
      "Weak Magnetism Neutrino Reactions", nuclib, opts, true, false, false,
      NeutrinoCorrectionMode::WeakMagnetism);
  ReactionNetwork narrowWmNetwork(nuclib, { &narrowWmReactionLib }, &helm,
      nullptr, opts);
  const double wmReferenceCenterMeV = 20.0;
  const double wmReferenceWidthMeV = 2.0;
  auto narrow20MeVHist = NarrowNeutrinoHistory(wmReferenceCenterMeV,
      wmReferenceWidthMeV, 1.e20);
  narrowWmNetwork.LoadNeutrinoHistory(narrow20MeVHist.MakeSharedPtr());
  std::vector<double> yFinalNarrowWm =
      narrowWmNetwork.Evolve(yInit, 0.0, 1.e10,
          FVT([&TGK] (double) {return TGK;}),
          FVT([&RHOGCC] (double) {return RHOGCC;}),
          "SkyNet_output_narrow_wm", 1.e-10).FinalY();

  const double electronMass = Constants::ElectronMassInMeV;
  const double eNueElectron = wmReferenceCenterMeV + deltaNP;
  const double eAntiNuePositron = wmReferenceCenterMeV - deltaNP;
  const double nueCaptureReference = eNueElectron
      * sqrt(eNueElectron * eNueElectron - electronMass * electronMass)
      * (1.0 + 1.01 * wmReferenceCenterMeV / 939.0);
  const double antiNueCaptureReference = eAntiNuePositron
      * sqrt(eAntiNuePositron * eAntiNuePositron
          - electronMass * electronMass)
      * (1.0 - 7.22 * wmReferenceCenterMeV / 939.0);
  const double approxYeFinNarrowWm = nueCaptureReference
      / (nueCaptureReference + antiNueCaptureReference);
  const double relativeErrNarrowWm = (approxYeFinNarrowWm
      - yFinalNarrowWm[1]) / yFinalNarrowWm[1];
  std::cout << yFinalNarrowWm[1] << " " << approxYeFinNarrowWm << " "
      << relativeErrNarrowWm << std::endl;

  if (fabs(relativeErrNarrowWm) > 2.e-3) {
    std::cerr << "WeakMagnetism narrow-spectrum equilibrium check failed: "
        << "corrected final proton abundance does not match the representative "
        << "20 MeV charged-current phase-space estimate with Horowitz "
        << "weak-magnetism factors." << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
