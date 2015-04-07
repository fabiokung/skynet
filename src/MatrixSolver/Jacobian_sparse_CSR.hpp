/// \file Jacobian_sparse_CSR.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_JACOBIAN_SPARSE_CSR_HPP_
#define SKYNET_MATRIXSOLVER_JACOBIAN_SPARSE_CSR_HPP_

#include <array>
#include <set>
#include <vector>

#include "MatrixSolver/Jacobian.hpp"
#include "Utilities/CodeError.hpp"

class Jacobian_sparse_CSR : public Jacobian {
public:
  Jacobian_sparse_CSR(const int n,
      const std::set<std::array<int, 2>>& nonZeroIdxs,
      MatrixSolver * const pMatrixSolver);

  void Reset() {
    for (unsigned int i = 0; i < mValues.size(); ++i)
      mValues[i] = 0.0;
  }

  double& operator()(const int i, const int j) {
    return mValues[FindFlatIndex(i, j)];
  }

  const double& operator()(const int i, const int j) const {
    return mValues[FindFlatIndex(i, j)];
  }

  void Describe(NetworkOutput * const pOutput) const;

  int NumRows() const {
    return mNumRows;
  }

  int NumNonZero() const {
    return (int)mValues.size();
  }

  const std::vector<int>& NumNonZeroPerRow() const {
    return mNumNonZeroPerRow;
  }

  const std::vector<int>& RowPtrs() const {
    return mRowPtrs;
  }

  const std::vector<int>& FORTRANRowPtrs() const {
    return mFORTRANRowPtrs;
  }

  const std::vector<int>& ColIdxs() const {
    return mColIdxs;
  }

  const std::vector<int>& FORTRANColIdxs() const {
    return mFORTRANColIdxs;
  }

  const std::vector<double>& Values() const {
    return mValues;
  }

private:
  int FindFlatIndex(const int i, const int j) const;

  int mNumRows;
  std::vector<int> mNumNonZeroPerRow;
  std::vector<int> mRowPtrs;
  std::vector<int> mFORTRANRowPtrs;
  std::vector<int> mColIdxs;
  std::vector<int> mFORTRANColIdxs;
  std::vector<double> mValues;
};

#endif // SKYNET_MATRIXSOLVER_JACOBIAN_SPARSE_CSR_HPP_
