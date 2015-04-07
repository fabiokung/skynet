/// \file Jacobian_lapack.hpp
/// \author jlippuner
/// \since Nov 29, 2014
///
/// \brief
///
///

#ifndef SKYNET_MATRIXSOLVER_JACOBIAN_LAPACK_HPP_
#define SKYNET_MATRIXSOLVER_JACOBIAN_LAPACK_HPP_

#include <cstring>

#include "MatrixSolver/Jacobian.hpp"

class Jacobian_lapack: public Jacobian {
public:
  Jacobian_lapack(const int n) :
      mN(n),
      mValues(n * n, 0.0) {}

  void Reset() {
    memset(mValues.data(), 0, mValues.size() * sizeof(double));
  }

  double& operator()(const int i, const int j) {
    return mValues[i + j * mN];
  }

  const double& operator()(const int i, const int j) const {
    return mValues[i + j * mN];
  }

  int N() const {
    return mN;
  }

  const std::vector<double>& Values() const {
    return mValues;
  }

  void Describe(NetworkOutput * const pOutput) const;

private:
  int mN;
  std::vector<double> mValues;
};

#endif // SKYNET_MATRIXSOLVER_JACOBIAN_LAPACK_HPP_
