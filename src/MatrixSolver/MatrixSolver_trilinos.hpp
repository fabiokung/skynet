/// \file MatrixSolver_trilinos.hpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_TRILINOS_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_TRILINOS_HPP_

#include <memory>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <Amesos.h>
#pragma GCC diagnostic pop

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

class MatrixSolver_trilinos : public MatrixSolver {
public:
  MatrixSolver_trilinos();

  void Analyze(const Jacobian * const pJacobian);

  std::vector<double> Solve(const Jacobian * const pJacobian,
      const std::vector<double>& rhs);

private:
  int mNRows;
  int mNNonZero;
  std::vector<int> mNumNonZeroPerRow;
  std::vector<int> mRowPtrs;
  std::vector<int> mColIdxs;

  std::unique_ptr<Epetra_CrsMatrix> mpMat;
  double * mpValues;
  std::unique_ptr<Epetra_Vector> mpRhs;
  std::unique_ptr<Epetra_Vector> mpResult;
  std::vector<int> mIndices;

  std::unique_ptr<Amesos_BaseSolver> mpSolver;
  Epetra_LinearProblem mProblem;
};


#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_TRILINOS_HPP_
