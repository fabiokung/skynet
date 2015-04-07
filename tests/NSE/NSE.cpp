/// \file NSE.cpp
/// \author jlippuner
/// \since Jul 30, 2014
///
/// \brief
///
///

#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Network/NSE.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "NuclearData/WebnucleoXmlReader.hpp"
#include "Utilities/Exp10.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  NSEOptions nseOpts;
  nseOpts.SetDoScreening(false);

  // n, p, he4
  {
    auto nuclidesLibrary =
        NuclideLibrary::CreateRestrictedLibrary(nuclib, { "n", "p", "he4" });
    NSE nse(nuclidesLibrary, &helm, nullptr, nseOpts);

//    for (int i = 0; i < nuclidesLibrary.NumNuclides(); ++i)
//      printf("%5s: BE = %.20E, m = %.20E\n", nuclidesLibrary.Names()[i].c_str(),
//          nuclidesLibrary.BindingEnergiesInMeV()[i],
//          nuclidesLibrary.MassesInMeV()[i]);

    auto nseResult = nse.CalcFromTemperatureAndDensity(5.0, 1.0E9, 0.1);
    auto nseY = nseResult.Y();
    std::vector<double> Yexpected = { 0.80000000000257665000,
        2.5766151590575760169E-12, 0.049999999998711693039 };

    double maxError = 0.0;
    printf("# final Y fraction, fractional error\n");
    for (unsigned int i = 0; i < nseY.size(); ++i) {
      if (std::isnan(nseY[i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }
      double error = fabs(nseY[i] - Yexpected[i]) / Yexpected[i];
      printf("# %5s: %10.4E %10.4E\n", nuclidesLibrary.Names()[i].c_str(),
          nseY[i], error);
      maxError = std::max(maxError, error);
    }

    if (maxError > 1.0E-10)
      return EXIT_FAILURE;
  }

  // NSE with X-Ray burst nuclides
  {
    std::vector<std::string> nuclideNames;
    std::vector<double> targetY;
    std::ifstream ifs("nse_xray_burst", std::ifstream::in);
    while (!ifs.eof()) {
      std::string name;
      double Y;
      ifs >> name >> Y;

      if (name.length() == 0)
        break;

      nuclideNames.push_back(name);
      targetY.push_back(Y);
    }
    ifs.close();

    auto nuclidesLibrary =
        NuclideLibrary::CreateRestrictedLibrary(nuclib, nuclideNames);
    NSE nse(nuclidesLibrary, &helm, nullptr, nseOpts);

    std::vector<double> nseY = nse.CalcFromTemperatureAndDensity(8.0, 1.0E6,
        0.25).Y();

    double maxError = 0.0;
    printf("# final Y fraction, absolute error\n");
    for (unsigned int i = 0; i < nseY.size(); ++i) {
      if (std::isnan(nseY[i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }

      double error = fabs(nseY[i] - targetY[i]);
      printf("# %5s: %10.4E %10.4E\n", nuclideNames[i].c_str(), nseY[i], error);
      maxError = std::max(maxError, error);
    }

    printf("max error = %g\n", maxError);
    if (maxError > 8.0E-4)
      return EXIT_FAILURE;
  }

  // NSE with all nuclides
  {
    std::vector<std::string> nuclideNames;
    std::vector<double> targetY;
    std::ifstream ifs("nse_all", std::ifstream::in);
    while (!ifs.eof()) {
      std::string name;
      double Y;
      ifs >> name >> Y;

      if (name.length() == 0)
        break;

      nuclideNames.push_back(name);
      targetY.push_back(Y);
    }
    ifs.close();

    auto nuclidesLibrary =
        NuclideLibrary::CreateRestrictedLibrary(nuclib, nuclideNames);
    NSE nse(nuclidesLibrary, &helm, nullptr, nseOpts);

    std::vector<double> nseY = nse.CalcFromTemperatureAndDensity(3.0, 1.0E8,
        0.15).Y();

    double maxError = 0.0;
    printf("# final Y fraction, absolute error\n");
    for (unsigned int i = 0; i < nseY.size(); ++i) {
      if (std::isnan(nseY[i])) {
        printf("got NaN\n");
        return EXIT_FAILURE;
      }

      double error = fabs(nseY[i] - targetY[i]);
      printf("# %5s: %10.4E %10.4E\n", nuclideNames[i].c_str(), nseY[i], error);
      maxError = std::max(maxError, error);
    }

    printf("max error = %g\n", maxError);
    if (maxError > 3.5E-5)
      return EXIT_FAILURE;
  }

  // check NSE for full parameter space
  {
    // temperature limits in GK
    double T9Low = 2.0;
    double T9High = 40.0;
    int numT9 = 6;

    // density limits in g / cm^3
    double rhoLow = 1.0;
    double rhoHigh = 1.0E12;
    int numRho = 4;

    // Ye limits
    double YeLow = 0.01;
    double YeHigh = 0.99;
    int numYe = 4;

    auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
        SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");
    NSE nse(nuclib, &helm, nullptr, nseOpts);

    NetworkOptions opts;
    REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "/data/reaclib",
        ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Strong reactions", nuclib, opts, false);

    printf("#%5s  %10s  %6s%6s  %10s  %10s  %10s\n", "T9", "rho", "Ye", "#it",
        "error", "inv error", "max yDot");
    for (int t = 0; t < numT9; ++t) {
      double T9 = T9Low + (T9High - T9Low) / (double)(numT9 - 1) * (double)t;

      for (int r = 0; r < numRho; ++r) {
        double rho = exp10(log10(rhoLow) + (log10(rhoHigh) - log10(rhoLow))
            / (double)(numRho - 1) * (double)r);

        for (int y = 0; y < numYe; ++y) {
          double Ye = YeLow
              + (YeHigh - YeLow) / (double)(numYe - 1) * (double)y;

          printf("%6.2f  %10.3E  %6.2f", T9, rho, Ye);
          auto res = nse.CalcFromTemperatureAndDensity(T9, rho, Ye);

          if (res.NumIterations() < 10) {
            // test getting NSE from entropy
            double entropy = helm.FromTAndRho(T9, rho, res.Y(), nuclib).S();
            double T9FromEntropy = nse.CalcFromEntropyAndDensity(entropy, rho,
                Ye).T9();

            printf("% 6i  %10.3E  %10.3E", res.NumIterations(), res.Error(),
                T9 - T9FromEntropy);

            if (fabs(T9 / T9FromEntropy - 1.0) > 1.0E-12) {
              printf("T9 / T9FromEntropy - 1.0 too large\n");
              return EXIT_FAILURE;
            }
          } else {
            printf("% 6i  %10.3E  %10s", res.NumIterations(), res.Error(),
                "skip inv");
          }

          // check NSE consistency, all strong reactions should give 0 y dots
          ThermodynamicState thermoState(T9, rho, 0.0, 0.0, 0.0, 0.0, 0.0);
          strongReactionLibrary.CalculateRates(thermoState,
              nuclib.PartitionFunctionsWithoutSpinTerm(T9,
                  opts.PartitionFunctionLog10Cap), opts.RateExpArgumentCap,
                  nullptr, nullptr);

          std::vector<double> yDot(res.Y().size(), 0.0);
          strongReactionLibrary.AddYdotContributions(res.Y(), &yDot);

          // find largest yDot
          double largestYDot = 0.0;
          for (unsigned int i = 0; i < yDot.size(); ++i) {
            if (fabs(yDot[i]) > largestYDot)
              largestYDot = fabs(yDot[i]);
          }

          printf("  %10.3E", largestYDot);

          if (largestYDot > 1.0E-200) {
            printf("  NSE consistency failed");
            return EXIT_FAILURE;
          }

          printf("\n");
        }
      }
    }
  }

  return EXIT_SUCCESS;
}
