/// \file LineSearch.cpp
/// \author jlippuner
/// \since Sep 20, 2014
///
/// \brief
///
///

#include "Utilities/LineSearch.hpp"

#include <cmath>
#include <stdexcept>

namespace { // unnamed so this can only be used in this file

static constexpr double Alpha = 1.0E-4;

}

std::vector<double> LineSearch::FindNextPoint(
    const std::vector<double>& xOld, const double gradFdotP,
    const std::vector<double>& p,
    const std::function<double(const std::vector<double>& x)>& calcF) {
  if (xOld.size() != p.size())
    throw std::invalid_argument("The vectors xOld and p must have the same "
        "size in LineSearch");

  if (gradFdotP >= 0.0)
    throw std::invalid_argument("positive gradFdotP in LineSearch>");

  double fOld = calcF(xOld);
  double lambdaMin = 0.1;
  double lambda = 1.0;
  double lambda2 = lambda;
  double f2 = fOld;

  std::vector<double> nextX(xOld.size());

  while (true) {
    for (unsigned int i = 0; i < xOld.size(); ++i)
      nextX[i] = xOld[i] + lambda * p[i];

    if (lambda < lambdaMin)
      return nextX;

    double f = calcF(nextX);

    if (f <= fOld + Alpha * lambda * gradFdotP)
      return nextX;

    double lambdaNew;
    if (lambda == 1.0) {
      // first backtrack
      lambdaNew = -gradFdotP / (2.0 * (f - fOld - gradFdotP));
    } else {
      // subsequent backtracks
      double rhs1 = f - fOld - lambda * gradFdotP;
      double rhs2 = f2 - fOld - lambda2 * gradFdotP;
      double divisor = 1.0 / (lambda - lambda2);
      double a = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2))
          * divisor;
      double b = (-lambda2 * rhs1 / (lambda * lambda)
          + lambda * rhs2 / (lambda2 * lambda2)) * divisor;

      if (a == 0.0) {
        lambdaNew = -gradFdotP / (2.0 * b);
      } else {
        double discriminant = b * b - 3.0 * a * gradFdotP;
        if (discriminant < 0.0)
          lambdaNew = 0.5 * lambda;
        else if (b <= 0.0)
          lambdaNew = (-b + sqrt(discriminant)) / (3.0 * a);
        else
          lambdaNew = -gradFdotP / (b + sqrt(discriminant));
      }

      lambdaNew = std::min(lambdaNew, 0.5 * lambda);
    }

    lambda2 = lambda;
    f2 = f;
    lambda = std::max(lambdaNew, 0.1 * lambda);
  }

  return nextX;
}
