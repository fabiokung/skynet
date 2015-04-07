/// \file SmallNetwork.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/TemperatureDensityHistory.hpp"
#include "Reactions/FFNReactionLibrary.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"
#include "Utilities/Profiler.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  // read nuclides and initial mass fractions
  std::vector<std::string> nuclideNames;
  std::vector<double> initialMassFraction;

  {
    std::ifstream ifs("initial_mass_fractions", std::ifstream::in);
    while (!ifs.eof()) {
      std::string name;
      double massFraction;

      ifs >> name >> massFraction;

      if (name.length() == 0)
        break;

      nuclideNames.push_back(name);
      initialMassFraction.push_back(massFraction);
    }

    ifs.close();
  }

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nuclideNames);

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  auto hist = TemperatureDensityHistory::CreateFromFile("temp_rho_vs_time");

  NetworkOptions opts;
  opts.MaxYChangePerStep = 1.0E-1;
  opts.SmallestYUsedForDtCalculation = 1.0E-6;
  opts.ConvergenceCriterion = NetworkConvergenceCriterionStruct::Mass;
  opts.MassDeviationThreshold = 1.0E-8;
  opts.SmallestYUsedForErrorCalculation = 1.0E-20;
  opts.MaxDtChangeMultiplier = 2.0;
  opts.EnableScreening = true;

  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, true);
  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, true);

  SkyNetScreening screen(nuclib);
  ReactionNetwork net(nuclib, { &weakReactionLibrary, &strongReactionLibrary},
      &helm, &screen, opts);

  std::vector<double> yInit(initialMassFraction.size());
  for (unsigned int i = 0; i < initialMassFraction.size(); ++i) {
    yInit[i] = initialMassFraction[i]
        / (double)net.GetNuclideLibrary().As().at(i);
  }

  auto output = net.Evolve(yInit, hist, "SkyNet_output");

  std::vector<double> yFinal = output.FinalY();

  std::vector<double> massFractionFinal = yFinal;
  for (unsigned int i = 0; i < massFractionFinal.size(); ++i) {
    massFractionFinal[i] *= (double)net.GetNuclideLibrary().As().at(i);
  }

//  FILE * fout = fopen("fin", "w");
//  for (unsigned int i = 0; i < massFractionFinal.size(); ++i)
//    fprintf(fout, "%.20E\n", massFractionFinal[i]);
//  fclose(fout);

  std::vector<double> expectedFinalMassfractions(massFractionFinal.size());
  {
    std::ifstream ifs("final_mass_fractions", std::ifstream::in);
    for (unsigned int i = 0; i < expectedFinalMassfractions.size(); ++i)
      ifs >> expectedFinalMassfractions[i];

    ifs.close();
  }

  double maxError = 0.0;
  for (unsigned int i = 0; i < massFractionFinal.size(); ++i) {
    if (std::isnan(massFractionFinal[i])) {
      printf("got NaN\n");
      return EXIT_FAILURE;
    }
    if (massFractionFinal[i] < opts.SmallestYUsedForErrorCalculation)
      continue;
    double error = fabs(massFractionFinal[i] - expectedFinalMassfractions[i])
        / expectedFinalMassfractions[i];
    maxError = std::max(maxError, error);
  }

  printf("max fractional error = %.10E\n", maxError);

  if (maxError < 1.0E-7) {
    return EXIT_SUCCESS;
  } else {
    //for (unsigned int i = 0; i<massFractionFinal.size();i++) {
    //  std::cout << nuclideNames[i] << " " << massFractionFinal[i] << " " << expectedFinalMassfractions[i] << "\n";
    //}
    return EXIT_FAILURE;
  }
}

