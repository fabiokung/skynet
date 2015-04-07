add_executable(MatrixSolverTest MatrixSolver.cpp)

add_test(MatrixSolver MatrixSolverTest)

target_link_libraries(MatrixSolverTest
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
