/// \file Jacobian_sparse_CSR.cpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#if defined(SKYNET_USE_PARDISO)                                                 \
    || defined(SKYNET_USE_MKL)                                                  \
    || defined(SKYNET_USE_CUDA)                                                 \
    || defined(SKYNET_USE_CUDA_55)                                              \
    || defined(SKYNET_USE_PETSC)                                                \
    || defined(SKYNET_USE_TRILINOS_KLU)                                         \
    || defined(SKYNET_USE_TRILINOS_SUPERLU)                                     \
    || defined(SKYNET_USE_TRILINOS_UMFPACK)                                     \
    || defined(SKYNET_USE_TRILINOS_LAPACK)

#include "MatrixSolver/MatrixSolver.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "MatrixSolver/Jacobian_sparse_CSR.hpp"

std::unique_ptr<Jacobian> Jacobian::Create(const int n,
    const std::set<std::array<int, 2>>& indicesOfNonZeroEntries,
    MatrixSolver * const pMatrixSolver) {
  return std::unique_ptr<Jacobian>(new Jacobian_sparse_CSR(n,
      indicesOfNonZeroEntries, pMatrixSolver));
}

Jacobian_sparse_CSR::Jacobian_sparse_CSR(const int n,
    const std::set<std::array<int, 2>>& nonZeroIdxs,
    MatrixSolver * const pMatrixSolver) :
    mNumRows(n),
    mNumNonZeroPerRow(n),
    mRowPtrs(n + 1),
    mFORTRANRowPtrs(n + 1),
    mColIdxs(0),
    mFORTRANColIdxs(0),
    mValues(0) {
  // add diagonal entries
  std::set<std::array<int, 2>> extendedNonZeroIdxs(nonZeroIdxs);
  for (int i = 0; i < n; ++i)
    extendedNonZeroIdxs.insert({{i, i}});

  mColIdxs = std::vector<int>(extendedNonZeroIdxs.size());
  mFORTRANColIdxs = std::vector<int>(extendedNonZeroIdxs.size());
  mValues = std::vector<double>(extendedNonZeroIdxs.size());

  // check that all indices are valid
  for (auto& idx : extendedNonZeroIdxs) {
    ASSERT(idx[0] < n, "Invalid row index");
    ASSERT(idx[1] < n, "Invalid row index");
  }

  // sort column indices by row
  std::vector<std::vector<int>> columnIdxsVsRow(mNumRows);
  for (auto& idxs : extendedNonZeroIdxs)
    columnIdxsVsRow[idxs[0]].push_back(idxs[1]);

  int currentFlatIdx = 0;

  // loop over rows
  for (int row = 0; row < mNumRows; ++row) {
    mRowPtrs[row] = currentFlatIdx;

    std::vector<int> columnIdxs = columnIdxsVsRow[row];
    // sort the column indices in this row
    std::sort(columnIdxs.begin(), columnIdxs.end());

    for (unsigned int i = 0; i < columnIdxs.size(); ++i)
      mColIdxs[currentFlatIdx++] = columnIdxs[i];
  }

  // there is an extra index in mStartingFlatIndexVsRow that points beyond
  // the last flat index
  mRowPtrs[mNumRows] = currentFlatIdx;
  ASSERT((unsigned int)currentFlatIdx == mValues.size(),
      "Unexpected last flat index");

  // count number of non-zero entries per row
  for (int row = 0; row < mNumRows; ++row)
    mNumNonZeroPerRow[row] = mRowPtrs[row + 1] - mRowPtrs[row];

  // make fortran indices (just add 1...)
  for (unsigned int i = 0; i < mRowPtrs.size(); ++i)
    mFORTRANRowPtrs[i] = mRowPtrs[i] + 1;
  for (unsigned int i = 0; i < mColIdxs.size(); ++i)
    mFORTRANColIdxs[i] = mColIdxs[i] + 1;

  // run analysis phase
  pMatrixSolver->Analyze(this);
}

void Jacobian_sparse_CSR::Describe(NetworkOutput * const pOutput) const {
  pOutput->Log("# Using sparse %i x %i Jacbian with %lu (%.2f%%) non-zero "
      "entries\n", mNumRows, mNumRows, mValues.size(),
      100.0 * (double)mValues.size() / (double)(mNumRows * mNumRows));
}

int Jacobian_sparse_CSR::FindFlatIndex(const int i, const int j) const {
  assert(i >= 0);
  assert(i < mNumRows);

  int startIdx = mRowPtrs[i];
  int endIdx = mRowPtrs[i + 1];

  for (int k = startIdx; k < endIdx; ++k) {
    int colIdx = mColIdxs[k];
    if (colIdx == j)
      return k;
  }

  throw std::out_of_range("Sparse Jacobian does not contain a non-zero entry "
      "(" + std::to_string(i) + ", " + std::to_string(j) + ")");
}

#endif // defined(SKYNET_USE_PARDISO)
       // || defined(SKYNET_USE_MKL)
       // || defined(SKYNET_USE_CUDA)
       // || defined(SKYNET_USE_CUDA_55)
       // || defined(SKYNET_USE_PETSC)
       // || defined(SKYNET_USE_TRILINOS_KLU)
       // || defined(SKYNET_USE_TRILINOS_SUPERLU)
       // || defined(SKYNET_USE_TRILINOS_UMFPACK)
       // || defined(SKYNET_USE_TRILINOS_LAPACK)
