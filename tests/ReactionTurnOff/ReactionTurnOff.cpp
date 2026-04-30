/// \file ReactionTurnOff.cpp
/// \author jlippuner
/// \since Jan 16, 2015
///
/// \brief
///
///

#include <cstdio>
#include <cmath>
#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "BuildInfo.hpp"
#include "DensityProfiles/ExpTMinus3.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Reactions/REACLIB.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/FlushStdout.hpp"
#include "Utilities/Profiler.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.MaxYChangePerStep = 0.1;
  opts.SmallestYUsedForDtCalculation = 1.0E-6;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.SmallestYUsedForErrorCalculation = 1.0E-20;
  opts.MaxDtChangeMultiplier = 2.0;
  opts.MinDt = 1.0E-8;
  opts.IsSelfHeating = true;

  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, false);
  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, false);
  REACLIBReactionLibrary symmetricFission(
      SkyNetRoot + "/data/netsu_panov_symmetric_0neut",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Symmetric neutron induced fission with 0 free neutrons", nuclib, opts,
      false);
  REACLIBReactionLibrary spontaneousFission(
      SkyNetRoot + "/data/netsu_sfis_Roberts2010rates",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Spontaneous fission", nuclib, opts, false);

  HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");

  ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary,
      &symmetricFission, &spontaneousFission }, &helmEOS, nullptr, opts);

  std::vector<double> initialY(net.GetNuclideLibrary().NumNuclides());
  {
    std::ifstream ifs("init_y", std::ifstream::in);
    for (unsigned int i = 0; i < initialY.size(); ++i)
      ifs >> initialY[i];

    ifs.close();
  }

  auto densityProfile = ExpTMinus3(3342929.205467419, 0.29 / 1000.0);

  auto output = net.EvolveSelfHeatingWithInitialTemperature(initialY,
      0.65581005515227098, 0.73215787940170385, 0.010318768360120788,
      &densityProfile, "SkyNet");

  std::vector<double> finalY = output.FinalY();

//  FILE * fout = fopen("fin", "w");
//  for (unsigned int i = 0; i < finalY.size(); ++i)
//    fprintf(fout, "%.20E\n", finalY[i]);
//  fclose(fout);

  // Read expected Y values token-by-token via strtod so that subnormal values
  // (which cause errno=ERANGE and set the stream's failbit) do not poison all
  // subsequent reads.
  std::vector<double> expectedY(finalY.size(), 0.0);
  {
    std::ifstream ifs("final_y", std::ifstream::in);
    std::string token;
    for (unsigned int i = 0; i < expectedY.size(); ++i) {
      if (!(ifs >> token)) break;
      errno = 0;
      char* end;
      double v = strtod(token.c_str(), &end);
      expectedY[i] = (errno == ERANGE) ? 0.0 : v;
    }
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
    // skip entries where the reference is also below the tracking threshold
    if (expectedY[i] < opts.SmallestYUsedForErrorCalculation)
      continue;
    double error = fabs(finalY[i] - expectedY[i])
        / expectedY[i];
    maxError = std::max(maxError, error);
  }

  printf("max fractional error = %.10E\n", maxError);

  if (maxError < 1.0E-2)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

