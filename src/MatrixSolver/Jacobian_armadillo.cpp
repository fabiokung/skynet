/// \file Jacobian_armadillo.cpp
/// \author jlippuner
/// \since Oct 29, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_ARMADILLO

#include "MatrixSolver/Jacobian_armadillo.hpp"

#include <array>
#include <memory>
#include <set>

std::unique_ptr<Jacobian> Jacobian::Create(const int n,
    const std::set<std::array<int, 2>>& /*indicesOfNonZeroEntries*/,
    MatrixSolver * const /*pMatrixSolver*/) {
  return std::unique_ptr<Jacobian>(new Jacobian_armadillo(n));
}

void Jacobian_armadillo::Describe(NetworkOutput * const pOutput) const {
  pOutput->Log("# Using dense %i x %i Jacbian\n", mMat.n_rows, mMat.n_cols);
}

#endif // SKYNET_USE_ARMADILLO
