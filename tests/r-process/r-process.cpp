/// \file r-process.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <cstdio>
#include <chrono>
#include <cmath>
#include <fstream>
#include <map>
#include <thread>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/NSE.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/TemperatureDensityHistory.hpp"
#include "Utilities/FlushStdout.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/Profiler.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.MaxYChangePerStep = 0.1;
  opts.SmallestYUsedForDtCalculation = 1.0E-4;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.SmallestYUsedForErrorCalculation = 1.0E-20;
  opts.MaxDtChangeMultiplier = 2.0;
  opts.MinDt = 1.0E-16;
  opts.IsSelfHeating = true;

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

  auto hist = TemperatureDensityHistory::CreateFromFile("temp_rho_vs_time");

  HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");
  SkyNetScreening screen(nuclib);

  ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary,
      &symmetricFission, &spontaneousFission }, &helmEOS, &screen, opts);

  //auto output = net.EvolveFromNSE(hist, "SkyNet_output");
  auto output = net.EvolveFromNSE(hist.StartTime(), hist.StartTime() + 1.E-10,
      &hist.TemperatureVsTime(), &hist.DensityVsTime(), hist.Ye(),
      "SkyNet_output");

  std::vector<double> finalY = output.FinalY();

//  FILE * fout = fopen("fin", "w");
//  for (unsigned int i = 0; i < finalY.size(); ++i)
//    fprintf(fout, "%.20E\n", finalY[i]);
//  fclose(fout);

  std::vector<double> expectedY(finalY.size());
  {
    std::ifstream ifs("final_abundances", std::ifstream::in);
    for (unsigned int i = 0; i < expectedY.size(); ++i)
      ifs >> expectedY[i];

    ifs.close();
  }

  double maxError = 0.0;
  for (unsigned int i = 0; i < finalY.size(); ++i) {
    if (std::isnan(finalY[i])) {
      printf("got NaN\n");
      return EXIT_FAILURE;
    }
    if (finalY[i] < opts.SmallestYUsedForErrorCalculation)
      continue;
    double error = fabs(finalY[i] - expectedY[i])
        / expectedY[i];
    maxError = std::max(maxError, error);
  }

  printf("max fractional error = %.10E\n", maxError);

  if (maxError < 1.0E-8)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

