/// \file r-process.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <cstdio>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

bool CheckMassFractions(const int numZones,
    const std::vector<std::vector<double>> finalCMassFractions,
    const double * const pFinalFortranMassFractions, const double YThreshold,
    const double errorThreshold) {
  bool ok = true;

  for (int z = 0; z < numZones; ++z) {
    double maxError = 0.0;
    for (unsigned int i = 0; i < finalCMassFractions[z].size(); ++i) {
      if (std::isnan(finalCMassFractions[z][i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }
      if (finalCMassFractions[z][i] < YThreshold)
        continue;

      double result = pFinalFortranMassFractions[i * numZones + z];
      double expected = finalCMassFractions[z][i];

      double error = fabs(result - expected) / expected;
      maxError = std::max(maxError, error);
    }

    printf("%i: max fractional error = %.10E\n", z, maxError);

    if (maxError > errorThreshold)
      ok = false;
  }

  return ok;
}

bool CheckHeating(const int numZones, const std::vector<double> cHeating,
    const double * const pFortranHeating, const double errorThreshold) {
  bool ok = true;

  for (int z = 0; z < numZones; ++z) {
    double expected = cHeating[z];
    double error = fabs(pFortranHeating[z] - expected) / expected;

    printf("%i: heating fractional error = %.10E\n", z, error);

    if (error > errorThreshold)
      ok = false;
  }

  return ok;
}

bool CheckEntropyChange(const int numZones, const std::vector<double> cEntropy,
    const double * const pFortranEntropy, const double errorThreshold) {
  bool ok = true;

  for (int z = 0; z < numZones; ++z) {
    double expected = cEntropy[z];
    double error = fabs(pFortranEntropy[z] - expected) / expected;

    printf("%i: entropy change fractional error = %.10E\n", z, error);

    if (error > errorThreshold)
      ok = false;
  }

  return ok;
}

extern "C"
int CheckResults(const int numZones, const int numIsotopes,
    const int numFortranTimes, const double * const pTotalHeating,
    const double * const pTotalEntropyChanges,
    const double * const pFinalMassFractions) {
  FloatingPointExceptions::Enable();

  if (numZones != 2)
    throw std::invalid_argument("Expected 2 zones");

  // read list of nuclides and initial compositions
  std::vector<std::string> nuclideNames;
  std::vector<std::vector<double>> initialY(numZones);
  std::ifstream finit("initial_compositions", std::ifstream::in);
  std::string line;
  while (std::getline(finit, line)) {
    if (line[0] == '#')
      continue;

    std::istringstream iss(line);
    std::string name;
    int A, Z;

    if (!(iss >> name >> A >> Z))
      throw std::runtime_error("Error reading initial_compositions");

    nuclideNames.push_back(name);

    for (int i = 0; i < numZones; ++i) {
      double X;
      if (!(iss >> X))
        throw std::runtime_error("Error reading initial_compositions");

      initialY[i].push_back(X / (double)A);
    }
  }
  finit.close();

  if ((int)nuclideNames.size() != numIsotopes)
    throw std::runtime_error("Found unexpected number of nuclides");

  // read temperature and density histories
  std::vector<double> times;
  std::vector<std::vector<double>> T9s(numZones);
  std::vector<std::vector<double>> rhos(numZones);
  std::ifstream fhist("temperature_density", std::ifstream::in);
  while (std::getline(fhist, line)) {
    if (line[0] == '#')
      continue;

    std::istringstream iss(line);
    double time;

    if (!(iss >> time))
      throw std::runtime_error("Error reading initial_compositions");

    times.push_back(time);

    for (int i = 0; i < numZones; ++i) {
      double T9;
      if (!(iss >> T9))
        throw std::runtime_error("Error reading initial_compositions");

      T9s[i].push_back(T9);
    }

    for (int i = 0; i < numZones; ++i) {
      double rho;
      if (!(iss >> rho))
        throw std::runtime_error("Error reading initial_compositions");

      rhos[i].push_back(rho);
    }
  }
  fhist.close();

  // make temperature and density functions, we make it piecewise constant to
  // emulate the operator splitting interface
  std::vector<double> extendedTimes;
  std::vector<std::vector<double>> extendedT9s(numZones);
  std::vector<std::vector<double>> extendedRhos(numZones);

  int numTimes = times.size();

  extendedTimes.resize(2 * numTimes - 2);
  extendedTimes[0] = times[0];
  extendedTimes[2 * numTimes - 3] = times[numTimes - 1];

  for (int i = 0; i < numZones; ++i) {
    extendedT9s[i].resize(2 * numTimes - 2);
    extendedT9s[i][0] = T9s[i][0];
    extendedT9s[i][2 * numTimes - 3] = T9s[i][numTimes - 1];

    extendedRhos[i].resize(2 * numTimes - 2);
    extendedRhos[i][0] = rhos[i][0];
    extendedRhos[i][2 * numTimes - 3] = rhos[i][numTimes - 1];
  }

  for (int i = 1; i < numTimes - 1; ++i) {
    extendedTimes[2*i - 1] = times[i] * 0.9999;
    extendedTimes[2*i] = times[i] * 1.0001;

    for (int j = 0; j < numZones; ++j) {
      extendedT9s[j][2*i - 1] = T9s[j][i-1];
      extendedT9s[j][2*i] = T9s[j][i];

      extendedRhos[j][2*i - 1] = rhos[j][i-1];
      extendedRhos[j][2*i] = rhos[j][i];
    }
  }

  std::vector<PiecewiseLinearFunction> T9VsTimes;
  std::vector<PiecewiseLinearFunction> rhoVsTimes;

  for (int i = 0; i < numZones; ++i) {
    T9VsTimes.push_back(PiecewiseLinearFunction(extendedTimes, extendedT9s[i],
        false));
    rhoVsTimes.push_back(PiecewiseLinearFunction(extendedTimes, extendedRhos[i],
        false));
  }

  auto nuclib =
      NuclideLibrary::CreateFromWebnucleoXML(
          SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nuclideNames);

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");
  SkyNetScreening screening(nuclib);

  NetworkOptions opts;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.MassDeviationThreshold = 1.0E-6;
  opts.MaxDt = 1.0E-8;
  opts.IsSelfHeating = false;
  opts.DisableStdoutOutput = true;

  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, true);
  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, true);
  REACLIBReactionLibrary symmetricFission(
      SkyNetRoot + "/data/netsu_panov_symmetric_0neut",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts,
      false);
  REACLIBReactionLibrary spontaneousFission(
      SkyNetRoot + "/data/netsu_sfis_Roberts2010rates",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Spontaneous fission", nuclib, opts, false);

  // first run normal evolution (not using operator splitting)
  std::vector<double> heatingNormalEvolution(numZones);
  std::vector<double> entropyChangeNormalEvolution(numZones);
  std::vector<std::vector<double>> finalMassFractionNormalEvolution(numZones);

  for (int i = 0; i < numZones; ++i) {
    ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary,
        &symmetricFission, &spontaneousFission }, &helm, &screening, opts);

    char buffer[1024];
    sprintf(buffer, "SkyNet_output_C_normal_%06d", i);
    auto output = net.Evolve(initialY[i], times[0], times[numTimes - 2],
        &(T9VsTimes[i]), &(rhoVsTimes[i]), std::string(buffer));

    std::vector<double> finalY = output.FinalY();
    finalMassFractionNormalEvolution[i] = finalY;

    for (unsigned int j = 0; j < output.A().size(); ++j)
      finalMassFractionNormalEvolution[i][j] *= output.A()[j];

    // get total heating
    heatingNormalEvolution[i] = 0.0;
    for (unsigned int j = 0; j < output.NumEntries(); ++j)
      heatingNormalEvolution[i] += output.DtVsTime()[j]
          * output.HeatingRateVsTime()[j];

    entropyChangeNormalEvolution[i] =
        output.EntropyVsTime()[output.NumEntries() - 1]
            - output.EntropyVsTime()[0];
  }

  // now run operator splitting evolution so that we compare apples to apples,
  // since Fortran ran operator splitting evolution
  std::vector<double> heatingOpSplit(numZones);
  std::vector<double> entropyChangesOpSplit(numZones);
  std::vector<std::vector<double>> finalMassFractionOpSplit(numZones);

  for (int i = 0; i < numZones; ++i) {
    ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary,
        &symmetricFission, &spontaneousFission }, &helm, &screening, opts);

    char buffer[1024];
    sprintf(buffer, "SkyNet_output_C_op_split_%06d", i);
    net.InitOperatorSplitEvolution(initialY[i], times[0], T9s[i][0], rhos[i][0],
        std::string(buffer));

    heatingOpSplit[i] = 0.0;
    entropyChangesOpSplit[i] = 0.0;
    for (int t = 0; t < numFortranTimes - 1; ++t) {
      auto res = net.OperatorSplitTakeStep(T9s[i][t], rhos[i][t], times[t+1]);
      heatingOpSplit[i] += res.Heating;
      entropyChangesOpSplit[i] += res.EntropyChange;
    }

    net.FinishOperatorSplitEvolution();

    auto output = net.Output();

    std::vector<double> finalY = output.FinalY();
    finalMassFractionOpSplit[i] = finalY;

    for (unsigned int j = 0; j < output.A().size(); ++j)
      finalMassFractionOpSplit[i][j] *= output.A()[j];
  }

  // compare Fortran results to operator splitting evolution, these should be
  // bit identical
  if (!CheckMassFractions(numZones, finalMassFractionOpSplit,
      pFinalMassFractions, opts.SmallestYUsedForErrorCalculation, 1.0E-14))
    return EXIT_FAILURE;

  if (!CheckHeating(numZones, heatingOpSplit, pTotalHeating, 1.0E-14))
    return EXIT_FAILURE;

  if (!CheckEntropyChange(numZones, entropyChangesOpSplit, pTotalEntropyChanges,
      1.0E-14))
    return EXIT_FAILURE;

  // compare Fortran results to normal evolution (expect larger deviations here)
  if (!CheckMassFractions(numZones, finalMassFractionNormalEvolution,
      pFinalMassFractions, opts.SmallestYUsedForErrorCalculation, 0.01))
    return EXIT_FAILURE;

  if (!CheckHeating(numZones, heatingNormalEvolution, pTotalHeating, 0.01))
    return EXIT_FAILURE;

  if (!CheckEntropyChange(numZones, entropyChangeNormalEvolution,
      pTotalEntropyChanges, 0.01))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

