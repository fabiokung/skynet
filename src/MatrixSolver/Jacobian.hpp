/// \file Jacobian.hpp
/// \author jlippuner
/// \since Jun 30, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_JACOBIAN_HPP_
#define SKYNET_MATRIXSOLVER_JACOBIAN_HPP_

#include <array>
#include <memory>
#include <set>

#include "MatrixSolver/MatrixSolver.hpp"
#include "Network/NetworkOutput.hpp"

class MatrixSolver;

class Jacobian {
public:
  virtual ~Jacobian() { }

  static std::unique_ptr<Jacobian> Create(const int n,
      const std::set<std::array<int, 2>>& indicesOfNonZeroEntries,
      MatrixSolver * const pMatrixSolver);

  virtual void Reset() =0;

  virtual double& operator()(const int i, const int j) =0;

  virtual const double& operator()(const int i, const int j) const =0;

  virtual void Describe(NetworkOutput * const pOutput) const =0;
};

#endif // SKYNET_MATRIXSOLVER_JACOBIAN_HPP_
