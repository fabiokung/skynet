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

  std::vector<double> target {
      4.005550255755044e-09,
      4.041710055765802e-18,
      1.558452308819913e-16,
      1.896768614981382e-19,
      1.142924132782445e-14,
      8.165069717350282e-09,
      4.529564077428202e-08,
      7.55760923685825e-08,
      5.790772326022599e-07,
      2.890779232498946e-09,
      3.230222429909035e-07,
      7.136600802226571e-05,
      0.01779010284965102 };

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

