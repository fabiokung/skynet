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
    ReactionNetwork network(nuclib, { &reactionLib }, &helm, opts);
    std::vector<double> yFinal =
        network.Evolve(yInit, 0.0, 1222.0, FVT([&TGK] (double) {return TGK;}),
            FVT([&RHOGCC] (double) {return RHOGCC;}), "SkyNet_output",
            NSEOptions(), 1.e-10).FinalY();
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
  return EXIT_SUCCESS;
}
