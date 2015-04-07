add_custom_command(
  OUTPUT NSE_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/NSE/nse_xray_burst .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/NSE/nse_all .
)

add_executable(NSE NSE.cpp NSE_test_input_files)

add_test(NSE NSE)

target_link_libraries(NSE
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
