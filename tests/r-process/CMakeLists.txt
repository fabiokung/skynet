add_custom_command(
  OUTPUT r-process_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/r-process/temp_rho_vs_time .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/r-process/final_abundances .
)

add_executable(r-process_test r-process.cpp r-process_test_input_files)

add_test(r-process r-process_test)

target_link_libraries(r-process_test
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
