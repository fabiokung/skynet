/// \file MatrixSolver_pardiso.cpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_PARDISO

#include "MatrixSolver/MatrixSolver_pardiso.hpp"

#include <type_traits>
#include <omp.h>

#include "MatrixSolver/Jacobian_sparse_CSR.hpp"
#include "Utilities/FloatingPointExceptions.hpp"
#include "Utilities/OutputMuter.hpp"

std::unique_ptr<MatrixSolver> MatrixSolver::Create() {
  return std::unique_ptr<MatrixSolver>(new MatrixSolver_pardiso());
}

namespace { // unnamed so these functions can only be used in this file

// using 1 thread seems to be fastest (at least for the X-ray burst test)
constexpr static int NumberOfOpenMPThreads = 1;

// PARDISO extern definitions

extern "C" void pardisoinit(
    void * const pt,
    const int * const mtype,
    const int * const solver,
    int * const iparam,
    double * const dparam,
    int * const error);

extern "C" void pardiso(
    void * const pt,
    const int * const maxfct,
    const int * const mnum,
    const int * const mtype,
    const int * const phase,
    const int * const n,
    const double * const a,
    const int * const ia,
    const int * const ja,
    const int * const perm,
    const int * const nrhs,
    int * const iparam,
    const int * const msglvl,
    const double * const b,
    double * const x,
    int * const error,
    double * const dparam);

/*
extern "C" void pardiso_chkmatrix(int * mtype, int * n, double * a, int * ia,
    int * ij, int * error);
extern "C" void pardiso_chkvec(int * n, int * nrhs, double * b, int * error);
extern "C" void pardiso_printstats(int * mtype, int * n, double * a, int * ia,
    int * ja, int * nrhs, double * b, int * error);
*/

void PardisoCheckError(const int error) {
  std::string prefix = "*** PARDISO ERROR *** ";
  switch (error) {
  case 0:
    return;
  case -1:
    throw std::runtime_error(prefix + "Input inconsistent.");
  case -2:
    throw std::runtime_error(prefix + "Not enough memory.");
  case -3:
    throw std::runtime_error(prefix + "Reordering problem.");
  case -4:
    throw std::runtime_error(prefix + "Zero pivot, numerical fact. or "
        "iterative refinement problem.");
  case -5:
    throw std::runtime_error(prefix + "Unclassified (internal) error.");
  case -6:
    throw std::runtime_error(prefix + "Preordering failed (matrix types 11, 13 "
        "only).");
  case -7:
    throw std::runtime_error(prefix + "Diagonal matrix problem.");
  case -8:
    throw std::runtime_error(prefix + "32-bit integer overflow problem.");
  case -10:
    throw std::runtime_error(prefix + "No license file pardiso.lic found.");
  case -11:
    throw std::runtime_error(prefix + "License is expired.");
  case -12:
    throw std::runtime_error(prefix + "Wrong username or hostname.");
  case -100:
    throw std::runtime_error(prefix + "Reached maximum number of "
        "Krylov-subspace iteration in iterative solver.");
  case -101:
    throw std::runtime_error(prefix  + "No sufficient convergence in "
        "Krylov-subspace iteration within 25 iterations.");
  case -102:
    throw std::runtime_error(prefix + "Error in Krylov-subspace iteration.");
  case -103:
    throw std::runtime_error(prefix + "Break-Down in Krylov-subspace "
        "iteration.");
  }

  throw std::runtime_error(prefix + "Unknown error");
}

} // namespace [unnamed]

MatrixSolver_pardiso::MatrixSolver_pardiso() :
    mN(0),
    mError(0) {
  OutputMuter::Mute();
  pardisoinit(mPt, &mMtype, &mSolver, mIparam, mDparam, &mError);
  OutputMuter::Unmute();
  try {
    PardisoCheckError(mError);
  } catch (const std::exception& ex) {
    CleanUp();
    throw;
  }

  // don't override with default values
  mIparam[0] = 1;

  // set number of cores
  mIparam[2] = NumberOfOpenMPThreads;
  omp_set_num_threads(NumberOfOpenMPThreads);
}

MatrixSolver_pardiso::~MatrixSolver_pardiso() {
  CleanUp();
}

void MatrixSolver_pardiso::Analyze(const Jacobian * const pJacobian) {
  const Jacobian_sparse_CSR * pJac = nullptr;

  try {
    pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);
    mN = pJac->NumRows();
    mRowPtrs = pJac->FORTRANRowPtrs();
    mColIdxs = pJac->FORTRANColIdxs();
  } catch (const std::exception& ex) {
    CleanUp();
    throw std::runtime_error("Error in MatrixSolver_pardiso::Analyze: "
        + std::string(ex.what()));
  }

  // run analysis
  RunPardiso(pJac->Values().data(), nullptr, 11, nullptr);
}

std::vector<double> MatrixSolver_pardiso::Solve(
    const Jacobian * const pJacobian, const std::vector<double>& rhs) {
  const Jacobian_sparse_CSR * pJac = nullptr;
  std::vector<double> result(rhs);

  try {
    pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);
    if (mN != pJac->NumRows())
      throw std::invalid_argument("Called matrix solver with Jacobian of "
          "different size.");
  } catch (const std::exception& ex) {
    CleanUp();
    throw std::runtime_error("Error in MatrixSolver_pardiso::Solve: "
        + std::string(ex.what()));
  }

  RunPardiso(pJac->Values().data(), rhs.data(), 23, result.data());

  return result;
}

void MatrixSolver_pardiso::RunPardiso(const double * const values,
    const double * const rhs, const int phase, double * const result) {
  try {
    // pardiso throws floating point exceptions, so disable them before calling
    // pardiso
    FloatingPointExceptions::Disable();
    pardiso(mPt, &mMaxfct, &mMnum, &mMtype, &phase, &mN, values,
        mRowPtrs.data(), mColIdxs.data(), mPerm, &mNrhs, mIparam, &mMsglvl,
        rhs, result, &mError, mDparam);
    FloatingPointExceptions::Enable();
    PardisoCheckError(mError);
  } catch (const std::exception& ex) {
    printf("Catch in RunPardiso: %s\n", ex.what());
    CleanUp();
    throw ex;
  }
}

void MatrixSolver_pardiso::CleanUp() {
  // release memory
  int phase = -1;
  pardiso(mPt, &mMaxfct, &mMnum, &mMtype, &phase, &mN, nullptr, mRowPtrs.data(),
      mColIdxs.data(), mPerm, &mNrhs, mIparam, &mMsglvl, nullptr, nullptr,
      &mError, mDparam);
  PardisoCheckError(mError);
}

#endif // SKYNET_USE_PARDISO
