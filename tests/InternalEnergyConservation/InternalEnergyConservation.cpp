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
#include "Reactions/ConstReactionLibrary.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Network/ReactionNetwork.hpp"
#include "DensityProfiles/ConstantFunction.hpp"
#include "Utilities/FloatingPointExceptions.hpp"
#include "Utilities/Profiler.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  {
    std::array<double, 24> partf;
    partf.fill(0.0);

    NuclideLibrary nuclib({ Nuclide(1, 2, 0.0, 0.0, partf, "A"),
        Nuclide(3, 6, 0.0, 0.0, partf, "B") }, "made up");

    NetworkOptions opts;
    opts.IsSelfHeating = true;
    opts.MaxDt = 1.0e-1;

    Reaction r0({"A"}, {"B"}, {3}, {1}, false, false, false, "dummy", nuclib);

    ConstReactionLibrary reactionLibrary(ReactionType::Strong,
        "Constant reactions", "made up", {r0}, {1.0}, nuclib, opts);

    HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");

    ReactionNetwork net(nuclib, { &reactionLibrary }, &helmEOS, nullptr, opts);

    double rho = 1.0E6;
    auto rho_vs_time = ConstantFunction(rho);

    net.EvolveSelfHeatingWithInitialTemperature({ 0.5, 0.0 }, 0,
        1.0E2, 1, &rho_vs_time, "toy_problem");

    // need to read HDF5 file because the output used in the network does
    // not keep all the Y's (to reduce memory usage)
    auto output = NetworkOutput::ReadFromFile("toy_problem.h5");

    // check energy is conserved
    double U0 = 0.0;
    double maxUError = 0.0;

    for (size_t i = 0; i < output.NumEntries(); ++i) {
      auto therm = helmEOS.FromTAndRho(output.TemperatureVsTime()[i], rho,
          output.YVsTime()[i], net.GetNuclideLibrary());
      double outU = therm.U();

      if (i == 0) U0 = outU;

      double thisUError = 2.0 * fabs(outU - U0) / fabs(outU + U0);
      maxUError = std::max(thisUError, maxUError);
    }

    if (maxUError > 3.0e-5) {
      printf("Internal energy not conserved in toy problem, max error "
          "= %.6e\n", maxUError);
      return EXIT_FAILURE;
    }
  }


  {
    auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
        SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml",
        { "n", "p", "d", "t", "he3", "he4", "he6", "he7", "he8", "li6", "li7",
          "li8", "be7", "be8", "b8" });

    NetworkOptions opts;
    opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    opts.MassDeviationThreshold = 1.0E-10;
    opts.MaxDt = 3.0E4;
    opts.DensityTimeScaleMultiplier = 0.1;
    opts.IsSelfHeating = true;
    opts.EnableScreening = true;

    REACLIBReactionLibrary weakReactionLibrary(SkyNetRoot + "/data/reaclib",
        ReactionType::Weak, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Weak reactions", nuclib, opts, true);
    REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
        ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Strong reactions", nuclib, opts, true);

    HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");
    SkyNetScreening screen(nuclib);

    ReactionNetwork net(nuclib, { &weakReactionLibrary,
        &strongReactionLibrary }, &helmEOS, &screen, opts);
    auto& lib = net.GetNuclideLibrary();

    std::vector<double> Y0(lib.NumNuclides(), 0.0);
    Y0[lib.NuclideIdsVsNames().at("p")] = 0.75;
    Y0[lib.NuclideIdsVsNames().at("he4")] = 0.25 / 4.0;

    double rho = 1.0E7;
    auto rho_vs_time = ConstantFunction(rho);

    net.EvolveSelfHeatingWithInitialTemperature(Y0, 0.0, 1.0E8, 3.0,
        &rho_vs_time, "h_burn");

    // need to read HDF5 file because the output used in the network does
    // not keep all the Y's (to reduce memory usage)
    auto output = NetworkOutput::ReadFromFile("h_burn.h5");

    // check energy is conserved
    double U0 = 0.0;
    double maxUError = 0.0;

    for (size_t i = 0; i < output.NumEntries(); ++i) {
      auto therm = helmEOS.FromTAndRho(output.TemperatureVsTime()[i], rho,
          output.YVsTime()[i], lib);
      double outU = therm.U();

      if (i == 0) U0 = outU;

      double thisUError = 2.0 * fabs(outU - U0) / fabs(outU + U0);
      maxUError = std::max(thisUError, maxUError);
    }

    if (maxUError > 5.0e-4) {
      printf("Internal energy not conserved in hydrogen burn, max error "
          "= %.6e\n", maxUError);
      return EXIT_FAILURE;
    }
  }

  {
    auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
        SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

    NetworkOptions opts;
    opts.ConvergenceCriterion = NetworkConvergenceCriterion::Mass;
    opts.MassDeviationThreshold = 1.0E-8;
    opts.MinDt = 1.0E-8;
    opts.MaxDt = 1.0E-2;
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

    double rho = 3.5E8;
    ConstantFunction rhoVsTime(rho);

    net.EvolveSelfHeatingFromNSEWithInitialEntropy(0.0025, 0.25, 7.7,
        &rhoVsTime, 0.1, "nse");

    // need to read HDF5 file because the output used in the network does
    // not keep all the Y's (to reduce memory usage)
    auto output = NetworkOutput::ReadFromFile("nse.h5");

    // check energy is conserved
    double U0 = 0.0;
    double maxUError = 0.0;

    FILE * f = fopen("U", "w");

    for (size_t i = 0; i < output.NumEntries(); ++i) {
      auto therm = helmEOS.FromTAndRho(output.TemperatureVsTime()[i],
          output.DensityVsTime()[i], output.YVsTime()[i],
          net.GetNuclideLibrary());
      double outU = therm.U();
      fprintf(f, "%20.10e  %20.10e\n", output.Times()[i], outU);

      if (i == 0)
        U0 = outU;

      double thisUError = 2.0 * fabs(outU - U0) / fabs(outU + U0);
      maxUError = std::max(thisUError, maxUError);
    }

    fclose(f);

    if (maxUError > 2.0e-05) {
      printf("Internal energy not conserved in NSE test, max error "
          "= %.6e\n", maxUError);
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
