/// \file MatrixSolver_mkl.cpp
/// \author jlippuner
/// \since Aug 2, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_MKL

#include "MatrixSolver/MatrixSolver_mkl.hpp"

#include <cmath>
#include <type_traits>

#include "MatrixSolver/Jacobian_sparse_CSR.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

std::unique_ptr<MatrixSolver> MatrixSolver::Create() {
  return std::unique_ptr<MatrixSolver>(new MatrixSolver_mkl());
}

#define DSS_CHK(value) {                                                        \
  MKL_INT status = value;                                                       \
  if (status != MKL_DSS_SUCCESS)                                                \
    throw std::runtime_error("** MKL DSS error at line "                        \
        + std::to_string(__LINE__) + " in file " + __FILE__);                   \
  }

MatrixSolver_mkl::MatrixSolver_mkl() :
    mNRows(0),
    mNCols(0),
    mNNonZero(0) {
  DSS_CHK(dss_create(mHandle, mOpt));
}

MatrixSolver_mkl::~MatrixSolver_mkl() {
  CleanUp();
}

void MatrixSolver_mkl::Analyze(const Jacobian * const pJacobian) {
  const Jacobian_sparse_CSR * pJac = nullptr;

  try {
    pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);
    mNRows = pJac->NumRows();
    mNCols = mNRows;
    mNNonZero = pJac->NumNonZero();
    mRowPtrs = pJac->FORTRANRowPtrs();
    mColIdxs = pJac->FORTRANColIdxs();
  } catch (const std::exception& ex) {
    CleanUp();
    throw std::runtime_error("Error in MatrixSolver_mkl::Analyze: "
        + std::string(ex.what()));
  }

  // define struction and reorder
  DSS_CHK(dss_define_structure(mHandle, mSym, mRowPtrs.data(), mNRows, mNCols,
      mColIdxs.data(), mNNonZero));
  DSS_CHK(dss_reorder(mHandle, mOpt, 0));
}

std::vector<double> MatrixSolver_mkl::Solve(
    const Jacobian * const pJacobian, const std::vector<double>& rhs) {
  const Jacobian_sparse_CSR * pJac = nullptr;
  std::vector<double> result(rhs);

  try {
    pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);
    if (mNRows != pJac->NumRows())
      throw std::invalid_argument("Called matrix solver with Jacobian of "
          "different size.");
  } catch (const std::exception& ex) {
    CleanUp();
    throw std::runtime_error("Error in MatrixSolver_mkl::Analyze: "
        + std::string(ex.what()));
  }

  DSS_CHK(dss_factor_real(mHandle, mType, pJac->Values().data()));
  DSS_CHK(dss_solve_real(mHandle, mOpt, rhs.data(), mNrhs, result.data()));

  return result;
}

void MatrixSolver_mkl::CleanUp() {
  // release memory
  DSS_CHK(dss_delete(mHandle, mOpt));
}

#endif // SKYNET_USE_MKL
