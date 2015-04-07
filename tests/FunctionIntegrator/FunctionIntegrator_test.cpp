/// \file FunctionIntegrator_test.cpp
/// \author jlippuner
/// \since Nov 26, 2014
///
/// \brief
///
///

#include <array>
#include <cmath>
#include <functional>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_fermi_dirac.h>

#include "Utilities/FunctionIntegrator.hpp"
#include "Utilities/StopWatch.hpp"

double WithSpecialFunctions(const double Q, const double mu) {
  // split result into 3 pieces
  double res = gsl_sf_fermi_dirac_int(4, mu) * gsl_sf_gamma(5.0);
  res -= 2.0 * Q * gsl_sf_fermi_dirac_int(3, mu) * gsl_sf_gamma(4.0);
  res += Q * Q * gsl_sf_fermi_dirac_int(2, mu) * gsl_sf_gamma(3.0);

  return res;
}

int main(int, char**) {
  // evaluate \int_0^\infty E^(E-Q)^2 / (1 + \exp(E-\mu))

  double minQ = 0.0;
  double maxQ = 10.0;
  int nQ = 101;

  double minMu = 0.0;
  double maxMu = 10.0;
  int nMu = 101;

  StopWatch swSpecialFunc;
  StopWatch swIntegrator;

  FunctionIntegrator integrator;

  printf("%12s  %12s  %12s  %12s  %12s\n", "Q", "mu", "spec func",
      "integrator", "error");

  double maxError = 0.0;

  for (int i = 0; i < nQ; ++i) {
    double Q = minQ + (maxQ - minQ) / (double)(nQ - 1) * (double)i;

    for (int j = 0; j < nMu; ++j) {
      double mu = minMu + (maxMu - minMu) / (double)(nMu - 1) * (double)j;

      swSpecialFunc.Start();
      double sf = WithSpecialFunctions(Q, mu);
      swSpecialFunc.Stop();

      swIntegrator.Start();
      auto in = integrator.Integrate([=] (const double E) {
            return E * E * (E - Q) * (E - Q) / (1.0 + exp(E - mu));
          }, 0.0, std::numeric_limits<double>::infinity());
      swIntegrator.Stop();

      double error = fabs(sf - in) / sf;
      printf("%12.5e  %12.5e  %12.5e  %12.5e  %12.5e\n", Q, mu, sf, in, error);

      maxError = std::max(maxError, error);
    }
  }

  printf("\n");
  printf("special function: %f ms\n",
      swSpecialFunc.GetElapsedTimeInMilliseconds());
  printf("      integrator: %f ms\n",
      swIntegrator.GetElapsedTimeInMilliseconds());
  printf("max error: %.5e\n", maxError);

  if (maxError < 1.0E-12)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}


