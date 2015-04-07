add_custom_command(
  OUTPUT XRayBurst_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/XRayBurst/initial_mass_fractions .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/XRayBurst/final_mass_fractions .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/XRayBurst/temp_rho_vs_time .
)

add_executable(XRayBurst XRayBurst.cpp XRayBurst_test_input_files)

add_test(XRayBurst XRayBurst)

target_link_libraries(XRayBurst
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
