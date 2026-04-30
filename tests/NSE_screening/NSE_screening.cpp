/// \file NSE_screening.cpp
/// \author jlippuner
/// \since Mar 21, 2017
///
/// \brief
///
///

#include <cmath>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "Utilities/Exp10.hpp"

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Network/NSE.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(SkyNetRoot
      + "data/webnucleo_nuc_v2.0.xml");

  HelmholtzEOS helm(SkyNetRoot + "data/helm_table.dat");
  SkyNetScreening screen(nuclib);

  NSEOptions nseOpts;
  nseOpts.SetDoScreening(true);

  // check NSE for full parameter space
  {
    // temperature limits in GK
    double T9Low = 5.0;
    double T9High = 10.0;
    int numT9 = 4;

    // density limits in g / cm^3
    double rhoLow = 1.0e8;
    double rhoHigh = 1.0e12;
    int numRho = 5;

    // Ye limits
    double YeLow = 0.01;
    double YeHigh = 0.99;
    int numYe = 5;

    NSE nse(nuclib, &helm, &screen, nseOpts);

    NetworkOptions opts;
    REACLIBReactionLibrary strongReactionLibrary(SkyNetRoot + "data/reaclib",
        ReactionType::Strong, true, LeptonMode::TreatAllAsDecayExceptLabelEC,
        "Strong reactions", nuclib, opts, nseOpts.DoScreening);

    int numFailed = 0;

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

          std::unique_ptr<NSEResult> res;
          try {
            res = std::unique_ptr<NSEResult>(
                new NSEResult(nse.CalcFromTemperatureAndDensity(T9, rho, Ye)));
          } catch (std::exception&) {
            printf(" failed\n");
            ++numFailed;
            // don't care for now FIXME
            //return EXIT_FAILURE;
            continue;
          }

          if (res->NumIterations() < 20) {
            // test getting NSE from entropy
            double entropy = helm.FromTAndRho(T9, rho, res->Y(), nuclib).S();
            double T9FromEntropy = nse.CalcFromEntropyAndDensity(entropy, rho,
                Ye).T9();

            printf("% 6i  %10.3E  %10.3E", res->NumIterations(), res->Error(),
                T9 - T9FromEntropy);

//            if (fabs(T9 / T9FromEntropy - 1.0) > 1.0E-6) {
//              printf("T9 / T9FromEntropy - 1.0 too large\n");
//              return EXIT_FAILURE;
//            }
          } else {
            printf("% 6i  %10.3E  %10s", res->NumIterations(), res->Error(),
                "skip inv");
          }

          // check NSE consistency, all strong reactions should give 0 y dots
          auto thermoState = helm.FromTAndRho(T9, rho, res->Y(), nuclib);
          auto scr = screen.Compute(thermoState, res->Y(), nuclib);
          auto mus = scr.Mus();


          strongReactionLibrary.CalculateRates(thermoState,
              nuclib.PartitionFunctionsWithoutSpinTerm(T9,
                  opts.PartitionFunctionLog10Cap), 250.0, //opts.RateExpArgumentCap,
                  &nuclib.Zs(), &mus);

          std::vector<double> yDot(res->Y().size(), 0.0);
          strongReactionLibrary.AddYdotContributions(res->Y(), &yDot);

          // find largest yDot
          double largestYDot = 0.0;
          for (unsigned int i = 0; i < yDot.size(); ++i) {
//            if (i < 10) printf("yDot[%u] = %10.3E\n", i, yDot[i]);
            if (fabs(yDot[i]) > largestYDot)
              largestYDot = fabs(yDot[i]);
          }

          printf("  %10.3E", largestYDot);

          if (largestYDot > 1.0E-20) {
            printf("  NSE consistency failed");
            return EXIT_FAILURE;
          }

          printf("\n");
//
//          for (unsigned int i = 0; i < yDot.size(); ++i) {
//            if (fabs(yDot[i]) > 1.0e-20)
//              printf("%6s: %12.6e\n", nuclib.Names()[i].c_str(), yDot[i]);
//          }
        }
      }
    }

    if (numFailed > 6)
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
