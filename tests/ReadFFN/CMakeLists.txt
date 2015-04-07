add_custom_command(
  OUTPUT ReadFFN_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/FFN/mesa_weak_rates.h5 .
)

add_executable(ReadFFN ReadFFN.cpp ReadFFN_test_input_files)

add_test(ReadFFN ReadFFN)

target_link_libraries(ReadFFN
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
