/// \file MatrixSolver_mkl.hpp
/// \author jlippuner
/// \since Aug 2, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_MKL_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_MKL_HPP_

#include <vector>

#include <mkl_dss.h>
#include <mkl_types.h>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

class MatrixSolver_mkl : public MatrixSolver {
public:
  MatrixSolver_mkl();

  ~MatrixSolver_mkl();

  void Analyze(const Jacobian * const pJacobian);

  std::vector<double> Solve(const Jacobian * const pJacobian,
      const std::vector<double>& rhs);

private:
  void CleanUp();

  MKL_INT mNRows;
  MKL_INT mNCols;
  MKL_INT mNNonZero;

  std::vector<MKL_INT> mRowPtrs;
  std::vector<MKL_INT> mColIdxs;

  _MKL_DSS_HANDLE_t mHandle;

  const MKL_INT mOpt = MKL_DSS_DEFAULTS | MKL_DSS_DEFAULTS;
  const MKL_INT mSym = MKL_DSS_NON_SYMMETRIC;
  const MKL_INT mType = MKL_DSS_INDEFINITE;
  const MKL_INT mNrhs = 1;
};

#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_MKL_HPP_
