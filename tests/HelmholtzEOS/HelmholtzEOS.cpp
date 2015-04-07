/// \file HelmholtzEOS.cpp
/// \author jlippuner
/// \since Aug 05, 2014
///
/// \brief
///
///

#include <cmath>
#include <cstdio>
#include <string>

#include "BuildInfo.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "NuclearData/WebnucleoXmlReader.hpp"
#include "Utilities/Exp10.hpp"
#include "Utilities/StopWatch.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  HelmholtzEOS helm(SkyNetRoot + "/data/helm_table.dat");

  std::vector<std::string> nuclideNames { "n", "p", "he4", "c12", "o16",
    "fe56" };

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml", nuclideNames);

  std::vector<double> Y{ 0.1, 0.2, 0.2 / 4.0, 0.2 / 12.0, 0.2 / 16.0,
    0.1 / 56.0 };

  // log10 of temperature (in Kelvin) limits
  int log10TLow = 4;
  int log10THigh = 13;

  // log10 of density (in g / cm^3) limits
  int log10RhoLow = -11;
  int log10RhoHigh = 12;

  StopWatch inversionStopWatch;
  double totalTime = 0.0;
  int numFailed = 0;
  double largestError = 0.0;

  printf("%14s  %14s  %14s  %14s  %14s\n", "T [GK]", "rho [cgs]", "S [kb/bry]",
      "frac error", "inv time [ms]");

  for (int log10T = log10TLow; log10T <= log10THigh; ++log10T) {
    for (int log10Rho = log10RhoLow; log10Rho <= log10RhoHigh; ++log10Rho) {
      double T9 = exp10((double)(log10T - 9) - 0.5);
      double rho = exp10((double)log10Rho);

      double entropy = helm.FromTAndRho(T9, rho, Y, nuclib).S();

      if (entropy < 1.0E-3) {
        printf("%14.6E  %14.6E  %14.6E  %s\n", T9, rho, entropy,
            "entropy too small, skip inversion");
      } else {
        inversionStopWatch.Start();
        double T9Inv = helm.FromSAndRho(entropy, rho, Y, nuclib).T9();
        inversionStopWatch.Stop();

        double error = helm.FromTAndRho(T9Inv, rho, Y, nuclib).S() / entropy - 1.0;

        if (fabs(error) > 1.0E-16)
          ++numFailed;

        largestError = std::max(largestError, fabs(error));

        double time = inversionStopWatch.GetElapsedTimeInMilliseconds();
        inversionStopWatch.Reset();
        totalTime += time;

        printf("%14.6E  %14.6E  %14.6E  %14.6E  %14.6F\n", T9, rho, entropy,
            error, time);
      }
    }
  }

  printf("total inversion time: %14.8F ms\n", totalTime);
  printf("number of failures: %i\n", numFailed);
  printf("largest error: %.4E\n", largestError);

  if ((numFailed > 30) || (largestError > 4.0E-8))
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
