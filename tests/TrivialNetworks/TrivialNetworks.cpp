/// \file TrivialNetwork.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <algorithm>
#include <cmath>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ConstReactionLibrary.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Utilities/FunctionVsTimeWrapper.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  {
    // dummy nuclides
    std::vector<Nuclide> nuclides = {
        Nuclide::CreateDummy(1, "A"),
        Nuclide::CreateDummy(1, "B")
    };
    NuclideLibrary nuclib(nuclides, "made up");

    // A -> B with r = 1
    Reaction r0({"A"}, {"B"}, {1}, {1}, true, false, false, "dummy", nuclib);
    // B -> A with r = 2
    Reaction r1({"B"}, {"A"}, {1}, {1}, true, false, false, "dummy", nuclib);

    NetworkOptions opts;
    opts.MaxDt = 0.1;

    ConstReactionLibrary reactionLibrary(ReactionType::Weak,
        "Constant reactions", "made up", {r0, r1}, {1.0, 2.0}, nuclib,
        opts);

    ReactionNetwork network(nuclib, { &reactionLibrary }, &helm, nullptr, opts);

    std::vector<double> yInit { 1.0, 0.0 };
    std::vector<double> yFinal =
        network.Evolve(yInit, 0.0, 10.0, FVT([] (double) { return 1.0; }),
            FVT([] (double) { return 1.0E9; }), "SkyNet_output").FinalY();

    std::vector<double> yAnalytic { 2.0 / 3.0, 1.0 / 3.0 };

    double maxError = 0.0;
    for (unsigned int i = 0; i < yFinal.size(); ++i) {
      if (std::isnan(yFinal[i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }
      double error = yFinal[i] - yAnalytic[i];
      maxError = std::max(maxError, fabs(error));
    }

    if (maxError > 1.0E-10) {
      printf("error %g is too big\n", maxError);
      return EXIT_FAILURE;
    }
  }

  {
    // dummy nuclides
    std::vector<Nuclide> nuclides = {
        Nuclide::CreateDummy(1, "A"),
        Nuclide::CreateDummy(1, "B"),
        Nuclide::CreateDummy(1, "C")
    };
    NuclideLibrary nuclib(nuclides, "made up");

    // cyclic network with 3 isotopes
    Reaction r0({"A"}, {"B"}, {1}, {1}, true, false, false, "dummy", nuclib);
    Reaction r1({"B"}, {"C"}, {1}, {1}, true, false, false, "dummy", nuclib);
    Reaction r2({"C"}, {"A"}, {1}, {1}, true, false, false, "dummy", nuclib);

    NetworkOptions opts;
    opts.MaxDt = 0.1;

    ConstReactionLibrary reactionLibrary(ReactionType::Weak,
        "Constant reactions", "made up", { r0, r1, r2 }, { 1.0, 1.5, 0.5 },
        nuclib, opts);

    ReactionNetwork network(nuclib, { &reactionLibrary }, &helm, nullptr, opts);

    std::vector<double> yInit { 1.0, 0.0, 0.0 };
    std::vector<double> yFinal =
        network.Evolve(yInit, 0.0, 10.0, FVT([] (double) { return 1.0; }),
            FVT([] (double) { return 1.0; }), "SkyNet_output").FinalY();

    std::vector<double> yAnalytic { 3.0 / 11.0, 2.0 / 11.0, 6.0 / 11.0 };

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

  // test that reactions that don't conserve mass are found
  {
    // dummy nuclides
    std::vector<Nuclide> nuclides = {
        Nuclide::CreateDummy(1, "A"),
        Nuclide::CreateDummy(2, "B")
    };
    NuclideLibrary nuclib(nuclides, "made up");

    // 2 A -> B with r = 1 (ok)
    Reaction r0( {"A"}, {"B"}, {2}, {1}, true, false, false, "dummy", nuclib);
    // B -> A with r = 2 (bad)
    Reaction r1( {"B"}, {"A"}, {1}, {1}, true, false, false, "dummy", nuclib);

    bool exceptionThrown = false;
    try {
      ConstReactionLibrary reactionLibrary0(ReactionType::Weak,
          "Constant reactions", "made up", { r0 }, { 1.0 },
          nuclib, NetworkOptions());
    } catch (...) {
      exceptionThrown = true;
    }

    if (exceptionThrown) {
      printf("falsely identified reaction that doesn't conserve mass\n");
      return EXIT_FAILURE;
    }

    exceptionThrown = false;
    try {
      ConstReactionLibrary reactionLibrary1(ReactionType::Weak,
          "Constant reactions", "made up", { r1 }, { 2.0 },
          nuclib, NetworkOptions());
    } catch (...) {
      exceptionThrown = true;
    }

    if (!exceptionThrown) {
      printf("didn't catch reaction that doesn't conserve mass\n");
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

