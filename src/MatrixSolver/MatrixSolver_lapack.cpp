/// \file MatrixSolver_lapack.cpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_LAPACK

#include "MatrixSolver/MatrixSolver_lapack.hpp"
#include "MatrixSolver/Jacobian_lapack.hpp"

extern "C" void dgesv_(
    const int * const N,      // number of equations
    const int * const NRHS,   // number of right-hand sides
    double * const A,         // values of matrix A (LDA-by-N, column-major)
    const int * const LDA,    // leading dimension of A (>= max(1,N))
    int * const IPIV,         // return array for pivots (N)
    double * const B,         // right-hand sides (LDB-by-NRHS, column-major)
    const int * const LDB,    // leading dimension of B (>= max(1,N))
    int * const INFO          // return 0 if successful
    );

std::unique_ptr<MatrixSolver> MatrixSolver::Create() {
  return std::unique_ptr<MatrixSolver>(new MatrixSolver_lapack());
}

void MatrixSolver_lapack::Analyze(const Jacobian * const /*pJacobian*/) {
  // do nothing
}

std::vector<double> MatrixSolver_lapack::Solve(
    const Jacobian * const pJacobian, const std::vector<double>& rhs) {
  auto pJac = dynamic_cast<const Jacobian_lapack * const>(pJacobian);

  // the LAPACK routine overwrites the input matrix and input rhs vector,
  // thus we need to copy them...

  int n = pJac->N();
  int nRhs = 1;
  std::vector<double> A(pJac->Values());
  std::vector<double> b(rhs);
  std::vector<int> pivots(n);

  int info;
  dgesv_(&n, &nRhs, A.data(), &n, pivots.data(), b.data(), &n, &info);

  if (info != 0)
    throw std::runtime_error("Error occurred in LAPACK DGESV routine, error "
        "code = " + std::to_string(info));

  return b;
}

#endif // SKYNET_USE_LAPACK

