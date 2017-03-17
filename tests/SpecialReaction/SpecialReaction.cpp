/// \file ReadFFN.cpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#include <vector>
#include <math.h>
#include <iostream>
#include "BuildInfo.hpp"
#include "Reactions/FFN.hpp"
#include "Utilities/FloatingPointExceptions.hpp"
#include "Reactions/SpecialReactionLibrary.hpp"
#include "Reactions/ReactionLibraryBase.hpp"
#include "Reactions/Reaction.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/EOS.hpp"
#include "Utilities/FunctionVsTimeWrapper.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/ReactionNetwork.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  std::vector<Nuclide> nuclides = {
      Nuclide::CreateDummy(1, "n"),
      Nuclide::CreateDummy(1, "p")
  };
  NuclideLibrary nuclib(nuclides, "Neutron Proton");

  NetworkOptions opts;
  opts.MaxDt = 0.001;

  const double rateConst = log(2.0)/1.0;

  auto reacNDecay =
      Reaction({"n"},{"p"},{1},{1}, true, false, true, "Neutron Decay", nuclib);

  RateFunc rateNDecay = [&rateConst](double /*time*/, ThermodynamicState state) {
    return rateConst*state.T9();
  };

  SpecialReactionLibrary specialLib(
      {std::pair<Reaction, RateFunc>(reacNDecay, rateNDecay)},
      ReactionType::Weak, "test", "Wikipedia", nuclib, opts);

  ReactionNetwork network(nuclib, { &specialLib }, &helm, opts);

  const double tfin = 2.0;
  const double beta = 0.2;
  std::vector<double> yInit { 1.0, 0.0 };
  std::vector<double> yFinal =
      network.Evolve(yInit, 0.0, tfin, FVT([beta] (double t) { return 1.0-beta*t; }),
          FVT([] (double) { return 1.0E9;}), "SkyNet_output").FinalY();

  // Compare answer to known value
  double numeric = yFinal[0];
  double analytic = exp(-rateConst * (tfin - 0.5*beta*tfin*tfin));
  std::cout << std::endl;
  std::cout << "Analytic : " << analytic << std::endl;
  std::cout << "Numeric : " << numeric << std::endl;
  if ( abs(analytic - numeric)/analytic < 0.01 ) return EXIT_SUCCESS;
  return EXIT_FAILURE;
}
