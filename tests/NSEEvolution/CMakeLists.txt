add_custom_command(
  OUTPUT NSEEvolution_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/NSEEvolution/particle .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/NSEEvolution/final_abundances .
)

add_executable(NSEEvolution NSEEvolution.cpp NSEEvolution_test_input_files)

add_test(NSEEvolution NSEEvolution)

target_link_libraries(NSEEvolution
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
