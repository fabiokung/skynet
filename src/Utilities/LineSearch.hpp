/// \file LineSearch.hpp
/// \author jlippuner
/// \since Sep 20, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_LINESEARCH_HPP_
#define SKYNET_UTILITIES_LINESEARCH_HPP_

#include <functional>
#include <vector>

/// Implements the line search and backtrack algorithm for Newton-Raphson
/// iterations described in section 9.7 of Numerical Recipes
class LineSearch {
public:
  static std::vector<double> FindNextPoint(const std::vector<double>& xOld,
      const double gradFdotP, const std::vector<double>& p,
      const std::function<double(const std::vector<double>& x)>& calcF);
};

#endif // SKYNET_UTILITIES_LINESEARCH_HPP_
