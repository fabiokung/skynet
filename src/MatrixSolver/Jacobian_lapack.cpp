/// \file Jacobian_lapack.cpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_LAPACK

#include "MatrixSolver/Jacobian_lapack.hpp"

#include <array>
#include <memory>
#include <set>

std::unique_ptr<Jacobian> Jacobian::Create(const int n,
    const std::set<std::array<int, 2>>& /*indicesOfNonZeroEntries*/,
    MatrixSolver * const /*pMatrixSolver*/) {
  return std::unique_ptr<Jacobian>(new Jacobian_lapack(n));
}

void Jacobian_lapack::Describe(NetworkOutput * const pOutput) const {
  pOutput->Log("# Using dense %i x %i Jacbian\n", mN, mN);
}

#endif // SKYNET_USE_LAPACK
