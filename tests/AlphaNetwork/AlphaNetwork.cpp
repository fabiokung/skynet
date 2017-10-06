/// \file SmallNetwork.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/FunctionVsTimeWrapper.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  std::vector<std::string> nuclideNames { "he4", "c12", "o16", "ne20", "mg24",
      "si28", "s32", "ar36", "ca40", "ti44", "cr48", "fe52", "ni56" };
  auto nuclideLibrary = NuclideLibrary::CreateFromWinv(
      SkyNetRoot + "/data/winvn_v2.0.dat", nuclideNames);

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  NetworkOptions opts;
  opts.MaxYChangePerStep = 5.0E-2;
  opts.SmallestYUsedForDtCalculation = 1.0E-10;
  opts.DeltaYByYThreshold = 1.0E-8;
  opts.SmallestYUsedForErrorCalculation = 1.0E-20;
  opts.MaxDtChangeMultiplier = 2.0;

  REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclideLibrary, opts, false);
  REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclideLibrary, opts, false);

  ReactionNetwork net(nuclideLibrary, { &weakReactionLibrary,
      &strongReactionLibrary }, &helm, nullptr, opts);

  std::vector<double> yInit { 0.0, 0.5 / 12.0, 0.5 / 16.0, 0.0, 0.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  std::vector<double> yFinal = net.Evolve(yInit, 1.0E-10, 1.0E12,
      FVT([] (double) {return 3.0;}), FVT([] (double) {return 1.0E9;}),
      "SkyNet_output").FinalY();

//  for (size_t i = 0; i < yFinal.size(); ++i) {
//    printf("%.15e\n", yFinal[i]);
//  }

  std::vector<double> target {
      4.227968124333936e-09,
      4.288114218916304e-18,
      9.120331634244609e-17,
      1.171603143449994e-19,
      7.418764420162448e-15,
      5.594320154167832e-09,
      3.275795161505568e-08,
      5.769221658804896e-08,
      4.665962942477050e-07,
      2.458619101403195e-09,
      2.899878322134674e-07,
      6.762531449578434e-05,
      1.779370524175675e-02 };

  double maxError = 0.0;
  printf("# final mass fraction, fractional error\n");
  for (unsigned int i = 0; i < yFinal.size(); ++i) {
    if (std::isnan(yFinal[i])) {
      printf("got NaN\n");
      return EXIT_FAILURE;
    }
    double error = fabs(yFinal[i] - target[i]) / target[i];
    printf("# %5s: %10.4E %10.4E\n", nuclideNames[i].c_str(), yFinal[i], error);
    maxError = std::max(maxError, error);
  }

  if (maxError < 5.0E-2)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

