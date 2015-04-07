/// \file MatrixSolver_armadillo.cpp
/// \author jlippuner
/// \since Oct 29, 2014
///
/// \brief
///
///

#ifdef SKYNET_USE_ARMADILLO

#include "MatrixSolver/MatrixSolver_armadillo.hpp"
#include "MatrixSolver/Jacobian_armadillo.hpp"

std::unique_ptr<MatrixSolver> MatrixSolver::Create() {
  return std::unique_ptr<MatrixSolver>(new MatrixSolver_armadillo());
}

void MatrixSolver_armadillo::Analyze(const Jacobian * const /*pJacobian*/) {
  // do nothing
}

std::vector<double> MatrixSolver_armadillo::Solve(
    const Jacobian * const pJacobian, const std::vector<double>& rhs) {;
  auto pJac = dynamic_cast<const Jacobian_armadillo * const>(pJacobian);
  arma::vec resultArma = arma::solve(pJac->Mat(), arma::vec(rhs));

  std::vector<double> result(resultArma.size());
  memcpy(result.data(), resultArma.mem, sizeof(double) * result.size());

  return  result;
}

#endif // SKYNET_USE_ARMADILLO
