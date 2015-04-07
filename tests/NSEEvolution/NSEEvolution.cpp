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
#include <sstream>
#include <thread>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "DensityProfiles/PowerLawContinuation.hpp"
#include "Utilities/FlushStdout.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/Profiler.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
  opts.MassDeviationThreshold = 1.0E-6;
  opts.MinDt = 1.0E-8;
  opts.DensityTimeScaleMultiplier = 0.1;
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

  HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");
  SkyNetScreening screen(nuclib);

  ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary,
      &symmetricFission, &spontaneousFission }, &helmEOS, &screen, opts);

  // read particle file
  double Ye = -1.0;
  bool foundYe = false;
  double initialEntropy = -1.0;
  bool foundInitialEntropy = false;
  std::vector<double> times;
  std::vector<double> densities;

  std::ifstream infile("particle", std::ifstream::in);
  std::string line;
  while (std::getline(infile, line)) {
    if (line[0] == '#') {
      // this is a comment, but it may contain Ye
      if (line.substr(0, 22) == "# Ye is constant, Ye =") {
        if (!foundYe) {
          std::string YeStr = line.substr(23);
          Ye = std::stod(YeStr);
          foundYe = true;
        } else {
          throw std::invalid_argument("Found second Ye value in particle file");
        }
      }

      continue;
    }

    // read values in line
    std::istringstream iss(line);
    double time, T9, rho, s;
    if (!(iss >> time >> T9 >> rho >> s))
      throw std::invalid_argument("Could not read line '" + line
          + "' in particle file");

    if (!foundInitialEntropy) {
      initialEntropy = s;
      foundInitialEntropy = true;
    }

    times.push_back(time);
    densities.push_back(rho);
  }

  infile.close();

  PiecewiseLinearFunction densityVsTime(times, densities, false);
  PowerLawContinuation densityProfile(densityVsTime, -3.0, 2.0);

  auto output = net.EvolveSelfHeatingFromNSEWithInitialEntropy(0.0025, 0.00335,
      initialEntropy, &densityProfile, Ye, "SkyNet_output");

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
    if (std::isnan (finalY[i])) {
      printf("got NaN\n");
      return EXIT_FAILURE;
    }
    if (finalY[i] < 1.0E-13)
      continue;

    double error = fabs(finalY[i] - expectedY[i]) / expectedY[i];
    //printf("%i: %.10E vs %.10E -> %.3E\n", i, finalY[i], expectedY[i], error);
    maxError = std::max(maxError, error);
  }

  printf("max fractional error = %.10E\n", maxError);

  if (maxError < 5.0E-5)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}
