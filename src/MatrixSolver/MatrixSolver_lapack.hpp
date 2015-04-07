/// \file MatrixSolver_lapack.hpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_LAPACK_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_LAPACK_HPP_

#include <vector>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

class MatrixSolver_lapack : public MatrixSolver {
public:
  void Analyze(const Jacobian * const /*pJacobian*/);

  std::vector<double> Solve(const Jacobian * const pJacobian,
      const std::vector<double>& rhs);
};

#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_LAPACK_HPP_
