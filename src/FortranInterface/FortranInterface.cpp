/// \file FortranInterface.cpp
/// \author jlippuner
/// \since Nov 1, 2014
///
/// \brief
///
///

#include <memory>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NSE.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"

// global unique pointers and variables that hold all necessary objects
NetworkOptions gNetworkOptions;
bool gDoScreening;
std::vector<std::unique_ptr<ReactionNetwork>> gpNetworks;
std::unique_ptr<NuclideLibrary> gpNuclib;
std::unique_ptr<EOS> gpEOS;
std::unique_ptr<Screening> gpScreening;
std::unique_ptr<NSE> gpNSE;

// global configuration functions, these must be called before SkyNetInit
extern "C"
void SkyNetSetConvergenceCriterion(const int convergenceCriterion) {
  switch (convergenceCriterion) {
  case 0:
    gNetworkOptions.ConvergenceCriterion =
        NetworkConvergenceCriterion::DeltaYByY;
    return;
  case 1:
    gNetworkOptions.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    return;
  case 2:
    gNetworkOptions.ConvergenceCriterion =
        NetworkConvergenceCriterion::BothDeltaYByYAndMass;
    return;
  }

  throw std::invalid_argument("Invalid number in "
      "SkyNetSetConvergenceCriterion. Specify one of the following: 1 for "
      "DeltaYByY, 2 for Mass, 3 for BothDeltaYByYAndMass.");
}

extern "C"
void SkyNetSetOptionDeltaYByYThreshold(const double deltaYByYThreshold) {
  gNetworkOptions.DeltaYByYThreshold = deltaYByYThreshold;
}

extern "C"
void SkyNetSetOptionMassDeviationThreshold(
    const double massDeviationThreshold) {
  gNetworkOptions.MassDeviationThreshold = massDeviationThreshold;
}

extern "C"
void SkyNetSetOptionMinDt(const double minDt) {
  gNetworkOptions.MinDt = minDt;
}

extern "C"
void SkyNetSetOptionMaxDt(const double maxDt) {
  gNetworkOptions.MaxDt = maxDt;
}

extern "C"
void SkyNetSetOptionDoScreening(const bool doScreening) {
  gDoScreening = doScreening;
}

// initialize nuclide library
extern "C"
void SkyNetInitNuclideLibrary(const int numIsotopes, const int * const pAs,
    const int * const pZs) {
  // assemble list of nuclides
  std::vector<std::string> isotopeNames(numIsotopes);
  for (int i = 0; i < numIsotopes; ++i)
    isotopeNames[i] = Nuclide::GetName(pZs[i], pAs[i]);

  //TODO add interface to give paths to these files
  auto tmpNuclib = NuclideLibrary::CreateFromWebnucleoXML(SkyNetRoot
      + "/data/webnucleo_nuc_v2.0.xml", isotopeNames);
  gpNuclib = std::unique_ptr<NuclideLibrary>(new NuclideLibrary(tmpNuclib));

  gpEOS = std::unique_ptr<HelmholtzEOS>(new HelmholtzEOS(SkyNetRoot
      + "/data/helm_table.dat", gNetworkOptions));
  gpScreening = std::unique_ptr<SkyNetScreening>(
      new SkyNetScreening(*gpNuclib.get()));

  // make NSE instance with the current nuclide library, we need to make a new
  // NSE instance after creating reaction networks, because they might change
  // the nuclide library
  gpNSE = std::unique_ptr<NSE>(new NSE(*gpNuclib.get(), gpEOS.get(),
      gpScreening.get()));
}

// initialize multiple copies of SkyNet
extern "C"
void SkyNetInit(const int numZones, const int numIsotopes,
    const double startTime, const double * const pInitialT9s,
    const double * const pInitialRhos, const int * const pAs,
    const int * const pZs, const double * const pInitialMassFractions) {
  // set options that are not settable by user
  gNetworkOptions.IsSelfHeating = false;
  gNetworkOptions.DisableStdoutOutput = true;
  //TODO option to disable h5 output

  SkyNetInitNuclideLibrary(numIsotopes, pAs, pZs);

  //TODO add interface to give paths to these files
  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", *gpNuclib.get(), gNetworkOptions, gDoScreening);
  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", *gpNuclib.get(), gNetworkOptions, gDoScreening);
  REACLIBReactionLibrary neutronFission(
      SkyNetRoot + "/data/netsu_panov_symmetric_0neut",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Neutron induced fission", *gpNuclib.get(), gNetworkOptions, false);
  REACLIBReactionLibrary spontaneousFission(
      SkyNetRoot + "/data/netsu_sfis_Roberts2010rates",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Spontaneous fission", *gpNuclib.get(), gNetworkOptions, false);

  // create and initialize networks
  gpNetworks.resize(numZones);
  for (int i = 0; i < numZones; ++i) {
    gpNetworks[i] = std::unique_ptr<ReactionNetwork>(
        new ReactionNetwork(*gpNuclib.get(), { &weakReactionLibrary,
            &strongReactionLibrary, &neutronFission, &spontaneousFission },
            gpEOS.get(), gDoScreening ? gpScreening.get() : nullptr,
                gNetworkOptions));

    // init network
    std::vector<double> initialY(numIsotopes);
    for (int j = 0; j < numIsotopes; ++j) {
      // the fortran array is arr(numZones, numIsotopes), which is stored in
      // column-major form
      initialY[j] = pInitialMassFractions[j * numZones + i] / pAs[j];
    }

    char buffer[1024];
    //TODO add interface to set output name
    sprintf(buffer, "SkyNet_output_Fortran_%06d", i);
    gpNetworks[i]->InitOperatorSplitEvolution(initialY, startTime,
        pInitialT9s[i], pInitialRhos[i], std::string(buffer));
  }

  // get new nuclide library
  if (numZones > 0) {
    gpNuclib = std::unique_ptr<NuclideLibrary>(
        new NuclideLibrary(gpNetworks[0]->GetNuclideLibrary()));

    gpNSE = std::unique_ptr<NSE>(new NSE(*gpNuclib.get(), gpEOS.get(),
        gpScreening.get()));
  }
}

// make SkyNets evolve at constant temperature and density until given end time,
// the specific heating rate and mass fractions at the end of the step are
// returned
extern "C"
void SkyNetTakeStep(const double * const pT9s, const double * const pRhos,
    const double endTime, double * const pSpecificHeating,
    double * const pEntropyChange, double * const pMassFractions) {
  for (unsigned int i = 0; i < gpNetworks.size(); ++i) {
    auto res = gpNetworks[i]->OperatorSplitTakeStep(pT9s[i], pRhos[i], endTime);
    pSpecificHeating[i] = res.Heating;
    pEntropyChange[i] = res.EntropyChange;

    for (unsigned int j = 0; j < gpNetworks[i]->Output().NumNuclides(); ++j) {
      // the fortran array is arr(numZones, numIsotopes), which is stored in
      // column-major form
      pMassFractions[j * gpNetworks.size() + i] =
          gpNetworks[i]->Output().FinalY()[j]
              * gpNetworks[i]->GetNuclideLibrary().As()[j];
    }
  }
}

// clean up SkyNets
extern "C"
void SkyNetFinish() {
  for (unsigned int i = 0; i < gpNetworks.size(); ++i)
    gpNetworks[i]->FinishOperatorSplitEvolution();

  // destroy networks
  gpNetworks.resize(0);

  // destroy nuclide library
  gpNuclib = nullptr;
}

/******************************************************************************
 **                                                                          **
 **  NSE INTERFACE                                                           **
 **                                                                          **
 **  MUST CALL SkyNetInitNuclideLibrary OR SkyNetInit FIRST                  **
 **                                                                          **
 ******************************************************************************/

void SkyNetNSEGetResult(const NSEResult& result, double * const pMassFractions,
    double * const pT9, double * const pRho) {
  for (int i = 0; i < gpNuclib->NumNuclides(); ++i)
    pMassFractions[i] = result.Y()[i] * gpNuclib->As()[i];

  if (pT9 != nullptr)
    *pT9 = result.T9();

  if (pRho != nullptr)
    *pRho = result.Rho();
}

extern "C"
void SkyNetNSECalcFromTemperatureAndDensity(const double T9, const double rho,
    const double Ye, double * const pMassFractions) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromTemperatureAndDensity(T9, rho, Ye),
      pMassFractions, nullptr, nullptr);
}

extern "C"
void SkyNetNSECalcFromTemperatureAndEntropy(const double T9,
    const double entropy, const double Ye, double * const pMassFractions,
    double * const pRho) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromTemperatureAndEntropy(T9, entropy, Ye), pMassFractions,
      nullptr, pRho);
}

extern "C"
void SkyNetNSECalcFromEntropyAndDensity(const double entropy, const double rho,
    const double Ye, double * const pMassFractions, double * const pT9) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromEntropyAndDensity(entropy, rho, Ye), pMassFractions,
      pT9, nullptr);
}

extern "C"
void SkyNetNSECalcFromEntropyAndDensityWithTemperatureGuess(
    const double entropy, const double rho, const double Ye,
    const double T9Guess, double * const pMassFractions, double * const pT9) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromEntropyAndDensityWithTemperatureGuess(entropy, rho, Ye,
          T9Guess), pMassFractions, pT9, nullptr);
}

extern "C"
void SkyNetNSECalcFromInternalEnergyAndDensity(const double internalEnergy,
    const double rho, const double Ye, double * const pMassFractions,
    double * const pT9) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromInternalEnergyAndDensity(internalEnergy, rho, Ye),
      pMassFractions, pT9, nullptr);
}

extern "C"
void SkyNetNSECalcFromInternalEnergyAndDensityWithTemperatureGuess(
    const double internalEnergy, const double rho, const double Ye,
    const double T9Guess, double * const pMassFractions, double * const pT9) {
  SkyNetNSEGetResult(
      gpNSE->CalcFromInternalEnergyAndDensityWithTemperatureGuess(
          internalEnergy, rho, Ye, T9Guess), pMassFractions, pT9, nullptr);
}

/******************************************************************************
 **                                                                          **
 **  EOS INTERFACE                                                           **
 **                                                                          **
 **  MUST CALL SkyNetInitNuclideLibrary OR SkyNetInit FIRST                  **
 **                                                                          **
 ******************************************************************************/

std::vector<double> SkyNetEOSMakeY(const double * const pMassFractions) {
  std::vector<double> Ys(gpNuclib->NumNuclides());
  for (int i = 0; i < gpNuclib->NumNuclides(); ++i)
    Ys[i] = pMassFractions[i] / gpNuclib->As()[i];

  return Ys;
}

extern "C"
void SkyNetEOSCalcFromTAndRho(const double T9, const double rho,
    const double * const pMassFractions, double * const pS, double * const pU,
    double * const pdSdT9) {
  auto res = gpEOS->FromTAndRho(T9, rho, SkyNetEOSMakeY(pMassFractions),
      *gpNuclib.get());

  *pS = res.S();
  *pU = res.U();
  *pdSdT9 = res.dSdT9();
}

extern "C"
void SkyNetEOSCalcFromSAndRho(const double entropy, const double rho,
    const double * const pMassFractions, double * const pT9, double * const pU,
    double * const pdSdT9) {
  auto res = gpEOS->FromSAndRho(entropy, rho, SkyNetEOSMakeY(pMassFractions),
      *gpNuclib.get());

  *pT9 = res.T9();
  *pU = res.U();
  *pdSdT9 = res.dSdT9();
}

extern "C"
void SkyNetEOSCalcFromSAndRhoWithTGuess(const double entropy, const double rho,
    const double T9Guess, const double * const pMassFractions,
    double * const pT9, double * const pU, double * const pdSdT9) {
  auto res = gpEOS->FromSAndRhoWithTGuess(entropy, rho,
      SkyNetEOSMakeY(pMassFractions), *gpNuclib.get(), T9Guess);

  *pT9 = res.T9();
  *pU = res.U();
  *pdSdT9 = res.dSdT9();
}

extern "C"
void SkyNetEOSCalcFromTAndS(const double T9, const double entropy,
    const double * const pMassFractions, double * const pRho, double * const pU,
    double * const pdSdT9) {
  auto res = gpEOS->FromTAndS(T9, entropy, SkyNetEOSMakeY(pMassFractions),
      *gpNuclib.get());

  *pRho = res.Rho();
  *pU = res.U();
  *pdSdT9 = res.dSdT9();
}

extern "C"
void SkyNetEOSCalcFromTAndSWithRhoGuess(const double T9, const double entropy,
    const double rhoGuess, const double * const pMassFractions,
    double * const pRho, double * const pU, double * const pdSdT9) {
  auto res = gpEOS->FromTAndSWithRhoGuess(T9, entropy,
      SkyNetEOSMakeY(pMassFractions), *gpNuclib.get(), rhoGuess);

  *pRho = res.Rho();
  *pU = res.U();
  *pdSdT9 = res.dSdT9();
}
