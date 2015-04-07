set(MatrixSolver_SOURCES
  Jacobian_armadillo.cpp
  Jacobian_lapack.cpp
  Jacobian_sparse_CSR.cpp
  MatrixSolver_armadillo.cpp
  MatrixSolver_lapack.cpp
  MatrixSolver_mkl.cpp
  MatrixSolver_pardiso.cpp
#  MatrixSolver_cuda.cpp
#  MatrixSolver_cuda_5.5.cpp
#  MatrixSolver_petsc.cpp
  MatrixSolver_trilinos.cpp
)

add_SkyNet_library(MatrixSolver "${MatrixSolver_SOURCES}")