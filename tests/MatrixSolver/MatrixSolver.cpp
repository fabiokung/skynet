/// \file MatrixSolver.cpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#include <array>
#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>

#include "MatrixSolver/Jacobian.hpp"
#include "MatrixSolver/Jacobian_sparse_CSR.hpp"
#include "MatrixSolver/MatrixSolver.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  std::set<std::array<int, 2>> nonZeroIndices{
    {{ 4, 1 }}, {{ 1, 1 }}, {{ 7, 2 }}, {{ 2, 2 }}, {{ 0, 2 }}, {{ 7, 6 }},
    {{ 0, 5 }}, {{ 1, 4 }}, {{ 5, 5 }}, {{ 7, 7 }}, {{ 5, 7 }}, {{ 2, 7 }},
    {{ 5, 2 }}, {{ 3, 3 }}, {{ 0, 0 }}, {{ 0, 6 }}, {{ 1, 2 }}, {{ 3, 6 }},
    {{ 6, 6 }}, {{ 6, 1 }}
  };

  // add some duplicates to the set
  nonZeroIndices.insert({{ 4, 1 }});
  nonZeroIndices.insert({{ 1, 4 }});
  nonZeroIndices.insert({{ 0, 0 }});
  nonZeroIndices.insert({{ 1, 2 }});
  nonZeroIndices.insert({{ 6, 1 }});

  auto solver = MatrixSolver::Create();
  auto jac = Jacobian::Create(8, nonZeroIndices, solver.get());
  jac->Reset();

  (*jac)(0, 6) = 7.0;
  (*jac)(0, 5) = 2.0;
  (*jac)(2, 7) = 5.0;
  (*jac)(1, 4) = 2.0;
  (*jac)(3, 6) = 9.0;
  (*jac)(0, 2) = 1.0;
  (*jac)(5, 7) = 5.0;
  (*jac)(1, 2) = 8.0;
  (*jac)(0, 0) = 7.0;
  (*jac)(1, 1) = -4.0;
  (*jac)(2, 2) = 1.0;
  (*jac)(3, 3) = 7.0;
  (*jac)(5, 5) = 3.0;
  (*jac)(6, 6) = 11.0;
  (*jac)(7, 7) = 5.0;
  (*jac)(7, 6) = 2.0;
  (*jac)(4, 1) = -4.0;
  (*jac)(5, 2) = 7.0;
  (*jac)(6, 1) = 17;
  (*jac)(7, 2) = -3.0;

#if defined(SKYNET_USE_PARDISO)                                                 \
    || defined(SKYNET_USE_MKL)                                                  \
    || defined(SKYNET_USE_CUDA)                                                 \
    || defined(SKYNET_USE_CUDA_55)                                              \
    || defined(SKYNET_USE_PETSC)                                                \
    || defined(SKYNET_USE_TRILINOS_KLU)                                         \
    || defined(SKYNET_USE_TRILINOS_SUPERLU)                                     \
    || defined(SKYNET_USE_TRILINOS_UMFPACK)                                     \
    || defined(SKYNET_USE_TRILINOS_LAPACK)

  Jacobian_sparse_CSR * pJacPardiso =
      dynamic_cast<Jacobian_sparse_CSR*>(jac.get());

  std::vector<int> rowPtrs { 0, 4, 7, 9, 11, 13, 16, 18, 21 };
  std::vector<int> colIdxs { 0, 2, 5, 6, 1, 2, 4, 2, 7, 3, 6, 1, 4, 2, 5, 7, 1,
    6, 2, 6, 7 };
  std::vector<double> values { 7.0, 1.0, 2.0, 7.0, -4.0, 8.0, 2.0, 1.0, 5.0,
      7.0, 9.0, -4.0, 0.0, 7.0, 3.0, 5.0, 17.0, 11.0, -3.0, 2.0, 5.0 };

  if (rowPtrs != pJacPardiso->RowPtrs())
    return EXIT_FAILURE;
  if (colIdxs != pJacPardiso->ColIdxs())
    return EXIT_FAILURE;
  if (values != pJacPardiso->Values())
    return EXIT_FAILURE;

  bool caughtException = false;
  try {
    (*jac)(7, 0) = 0.0;
  } catch (std::out_of_range&) {
    caughtException = true;
  }
  if (!caughtException) {
    printf("Failed to catch out of range exception in sparse matrix\n");
    return EXIT_FAILURE;
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

  std::vector<double> b { 8.0, 1.0, 6.0, -7.0, 4.0, -5.0, -7.0, 6.0 };
  std::vector<double> xAnalytic = { 31.0 / 21.0, -1.0, 5.0 / 11.0, -167.0
      / 77.0, -73.0 / 22.0, -151.0 / 33.0, 10.0 / 11.0, 61.0 / 55.0 };

  std::vector<double> x = solver->Solve(jac.get(), b);

  double maxError = 0.0;
  for (unsigned int i = 0; i < xAnalytic.size(); ++i) {
    if (std::isnan (x[i])) {
      printf("got NaN\n");
      return EXIT_FAILURE;
    }
    double error = fabs(x[i] - xAnalytic[i]) / xAnalytic[i];
    maxError = std::max(maxError, error);
  }

  printf("Max error = %.10e\n", maxError);

  if (maxError < 1.0E-10)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}


