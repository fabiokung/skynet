/// \file MatrixSolver_armadillo.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_ARMADILLO_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_ARMADILLO_HPP_

#include <armadillo>
#include <type_traits>
#include <vector>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

class MatrixSolver_armadillo : public MatrixSolver {
public:
  void Analyze(const Jacobian * const /*pJacobian*/);

  std::vector<double> Solve(const Jacobian * const pJacobian,
      const std::vector<double>& rhs);
};

#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_ARMADILLO_HPP_
