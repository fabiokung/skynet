add_custom_command(
  OUTPUT FortranOperatorSplitting_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/FortranOperatorSplitting/initial_compositions .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/FortranOperatorSplitting/temperature_density .
)

add_library(CheckResults CheckResults.cpp)

add_executable(FortranOperatorSplitting
  FortranOperatorSplitting.f90
  FortranOperatorSplitting_test_input_files
  ${PROJECT_SOURCE_DIR}/src/FortranInterface/SkyNet.f90
)

set_target_properties(FortranOperatorSplitting
  PROPERTIES LINKER_LANGUAGE Fortran
)

add_dependencies(FortranOperatorSplitting CheckResults)

add_test(FortranOperatorSplitting FortranOperatorSplitting)

target_link_libraries(FortranOperatorSplitting
  CheckResults
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
