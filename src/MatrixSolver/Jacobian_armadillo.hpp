/// \file Jacobian_armadillo.hpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_JACOBIAN_ARMADILLO_HPP_
#define SKYNET_MATRIXSOLVER_JACOBIAN_ARMADILLO_HPP_

#include <armadillo>

#include "MatrixSolver/Jacobian.hpp"

class Jacobian_armadillo: public Jacobian {
public:
  Jacobian_armadillo(const int n) :
      mMat(n, n, arma::fill::zeros) { }

  void Reset() {
    mMat.fill(0.0);
  }

  double& operator()(const int i, const int j) {
    return mMat(i, j);
  }

  const double& operator()(const int i, const int j) const {
    return mMat(i, j);
  }

  const arma::mat& Mat() const {
    return mMat;
  }

  void Describe(NetworkOutput * const pOutput) const;

private:
  arma::mat mMat;
};

#endif // SKYNET_MATRIXSOLVER_JACOBIAN_ARMADILLO_HPP_
