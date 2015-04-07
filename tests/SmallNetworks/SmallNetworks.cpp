/// \file SmallNetwork.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <algorithm>
#include <cmath>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ConstReactionLibrary.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/NetworkOptions.hpp"
#include "Utilities/FunctionVsTimeWrapper.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  {
    std::vector<Nuclide> nuclides = {
        Nuclide::CreateDummy(1, "A"),
        Nuclide::CreateDummy(1, "B"),
        Nuclide::CreateDummy(1, "C"),
        Nuclide::CreateDummy(2, "D") };
    NuclideLibrary nuclib(nuclides, "made up");

    // A + B <-> 2 C
    Reaction r0({"A", "B"}, {"C"}, {1, 1}, {2}, false, false, false, "dummy",
        nuclib);
    Reaction r0inv({"C"}, {"A", "B"}, {2}, {1, 1}, false, false, false,
        "dummy", nuclib);

    // 2 B + C <-> A + D
    Reaction r1({"B", "C"}, {"A", "D"}, {2, 1}, {1, 1}, false, false, false,
        "dummy", nuclib);
    Reaction r1inv({"A", "D"}, {"B", "C"}, {1, 1}, {2, 1}, false, false, false,
        "dummy", nuclib);

    // 3 C <-> B + D
    Reaction r2({"C"}, {"B", "D"}, {3}, {1, 1}, false, false, false, "dummy",
        nuclib);
    Reaction r2inv({"B", "D"}, {"C"}, {1, 1}, {3}, false, false, false,
        "dummy", nuclib);

    NetworkOptions opts;
    opts.MaxDt = 500.0;

    ConstReactionLibrary reactionLibrary(ReactionType::Strong,
        "Constant reactions", "made up", { r0, r0inv, r1, r1inv, r2, r2inv },
        { 1.0, 5.0, 5.0, 1.0, 0.1, 10.0 }, nuclib, opts);

    ReactionNetwork network(nuclib, { &reactionLibrary }, &helm, nullptr, opts);

    std::vector<double> yInit { 0.5, 0.1, 0.1, 0.15 };
    std::vector<double> yFinal = network.Evolve(yInit, 0.0, 15000.0,
        FVT([] (double) { return 1.0; }), FVT([] (double) { return 1.0; }),
        "SkyNet_output").FinalY();

    std::vector<double> yAnalytic { 0.91137, 0.0148826, 0.0736575,
      0.0000447529 };

    double maxError = 0.0;
    for (unsigned int i = 0; i < yFinal.size(); ++i) {
      if (std::isnan(yFinal[i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }
      double error = yFinal[i] - yAnalytic[i];
      maxError = std::max(maxError, fabs(error));
    }

    if (maxError > 1.0E-6) {
      printf("error %g is too big\n", maxError);
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}




