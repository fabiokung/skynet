/// \file MatrixSolver.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_HPP_

#include <memory>
#include <vector>

#include "MatrixSolver/Jacobian.hpp"

class Jacobian;

class MatrixSolver {
public:
  static std::unique_ptr<MatrixSolver> Create();

  virtual ~MatrixSolver() { }

  virtual void Analyze(const Jacobian * const pJacobian) =0;

  virtual std::vector<double> Solve(const Jacobian * const jacobian,
      const std::vector<double>& rhs) =0;
};

#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_HPP_
