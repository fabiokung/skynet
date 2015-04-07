/// \file MatrixSolver_pardiso.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_MATRIXSOLVER_PARDISO_HPP_
#define SKYNET_MATRIXSOLVER_MATRIXSOLVER_PARDISO_HPP_

#include <vector>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/MatrixSolver.hpp"

class MatrixSolver_pardiso : public MatrixSolver {
public:
  MatrixSolver_pardiso();

  ~MatrixSolver_pardiso();

  void Analyze(const Jacobian * const pJacobian);

  std::vector<double> Solve(const Jacobian * const pJacobian,
      const std::vector<double>& rhs);

private:
  void RunPardiso(const double * const values, const double * const rhs,
      const int phase, double * const result);

  void CleanUp();

  int mN;
  std::vector<int> mRowPtrs;
  std::vector<int> mColIdxs;

  void * mPt[64];
  const int mMaxfct = 1;
  const int mMnum = 1;
  const int mMtype = 11;
  const int mSolver = 0;
  const int * mPerm = nullptr;
  const int mNrhs = 1;
  int mIparam[64];
  const int mMsglvl = 0;
  double mDparam[64];
  int mError;
};

#endif // SKYNET_MATRIXSOLVER_MATRIXSOLVER_PARDISO_HPP_
