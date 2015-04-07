/// \file MatrixSolver_trilinos.cpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#if defined(SKYNET_USE_TRILINOS_KLU)                                            \
    || defined(SKYNET_USE_TRILINOS_SUPERLU)                                     \
    || defined(SKYNET_USE_TRILINOS_UMFPACK)                                     \
    || defined(SKYNET_USE_TRILINOS_LAPACK)

#include <Amesos_ConfigDefs.h>
#include <Epetra_SerialComm.h>
#include <Epetra_LocalMap.h>

#include "MatrixSolver/MatrixSolver_trilinos.hpp"
#include "MatrixSolver/Jacobian_sparse_CSR.hpp"

std::unique_ptr<MatrixSolver> MatrixSolver::Create() {
  return std::unique_ptr<MatrixSolver>(new MatrixSolver_trilinos());
}

#define AMESOS_CHK(value) {                                                     \
  if (value != 0)                                                               \
    throw std::runtime_error("** AMESOS error " + std::to_string(value)         \
        + " at line " + std::to_string(__LINE__) + " in file " + __FILE__);     \
  }

#define EPETRA_CHK(value) {                                                     \
  if (value != 0)                                                               \
    throw std::runtime_error("** EPETRA error " + std::to_string(value)         \
        + " at line " + std::to_string(__LINE__) + " in file " + __FILE__);     \
  }

MatrixSolver_trilinos::MatrixSolver_trilinos():
  mNRows(0) {}

void MatrixSolver_trilinos::Analyze(const Jacobian * const pJacobian) {
  const Jacobian_sparse_CSR * pJac = nullptr;

  pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);

  mNRows = pJac->NumRows();
  mNNonZero = pJac->NumNonZero();
  mNumNonZeroPerRow = pJac->NumNonZeroPerRow();
  mRowPtrs = pJac->RowPtrs();
  mColIdxs = pJac->ColIdxs();

  Epetra_LocalMap rowMap(mNRows, 0, Epetra_SerialComm());
  mpMat = std::unique_ptr<Epetra_CrsMatrix>(new Epetra_CrsMatrix(
      Epetra_DataAccess::Copy, rowMap, mNumNonZeroPerRow.data(), true));

  // fill with 1's so values don't get lost
  for (int i = 0; i < mNRows; ++i) {
    int nonZero = mNumNonZeroPerRow[i];
    std::vector<double> values(nonZero, 1.0);
    EPETRA_CHK(mpMat->InsertGlobalValues(i, mNumNonZeroPerRow[i],
        values.data(), mColIdxs.data() + mRowPtrs[i]));
  }

  EPETRA_CHK(mpMat->FillComplete());
  //EPETRA_CHK(mpMat->OptimizeStorage());
  //EPETRA_CHK(mpMat->MakeDataContiguous());

  // set column indices
  int * pRowPtrs;
  int * pColIdxs;
  if (mpMat->ExtractCrsDataPointers(pRowPtrs, pColIdxs, mpValues) != 0)
    throw std::runtime_error("ExtractCrsDataPointers failed");
  memcpy(pColIdxs, pJac->ColIdxs().data(), mNNonZero * sizeof(int));

  // make vectors
  mpRhs = std::unique_ptr<Epetra_Vector>(new Epetra_Vector(rowMap));
  mpResult = std::unique_ptr<Epetra_Vector>(new Epetra_Vector(rowMap));

  mIndices.resize(mNRows);
  for (int i = 0; i < mNRows; ++i)
    mIndices[i] = i;

  mProblem.SetOperator(mpMat.get());
  mProblem.SetLHS(mpResult.get());
  mProblem.SetRHS(mpRhs.get());

  Amesos factory;
  std::string solverName;

#if defined(SKYNET_USE_TRILINOS_KLU)
  solverName = "Amesos_Klu";
#elif defined(SKYNET_USE_TRILINOS_SUPERLU)
  solverName = "Amesos_Superlu";
#elif defined(SKYNET_USE_TRILINOS_UMFPACK)
  solverName = "Amesos_Umfpack";
#elif defined(SKYNET_USE_TRILINOS_LAPACK)
  solverName = "Amesos_Lapack";
#else
  solverName = "";
#endif

  mpSolver = std::unique_ptr<Amesos_BaseSolver>(factory.Create(solverName,
      mProblem));

  AMESOS_CHK(mpSolver->SymbolicFactorization());
}

std::vector<double> MatrixSolver_trilinos::Solve(
    const Jacobian * const pJacobian, const std::vector<double>& rhs) {
  const Jacobian_sparse_CSR * pJac = nullptr;
  std::vector<double> result(rhs);

  pJac = dynamic_cast<const Jacobian_sparse_CSR *>(pJacobian);
  if (mNRows != pJac->NumRows())
    throw std::invalid_argument("Called matrix solver with Jacobian of "
        "different size.");

  // put values into matrix
  memcpy(mpValues, pJac->Values().data(), mNNonZero * sizeof(double));

  EPETRA_CHK(mpRhs->ReplaceGlobalValues(mNRows, rhs.data(), mIndices.data()));

  AMESOS_CHK(mpSolver->NumericFactorization());
  AMESOS_CHK(mpSolver->Solve());

  EPETRA_CHK(mpResult->ExtractCopy(result.data()));
  return result;
}

#endif // defined(SKYNET_USE_TRILINOS_KLU)
       // || defined(SKYNET_USE_TRILINOS_SUPERLU)
       // || defined(SKYNET_USE_TRILINOS_UMFPACK)
       // || defined(SKYNET_USE_TRILINOS_LAPACK)
