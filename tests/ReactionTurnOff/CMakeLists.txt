add_custom_command(
  OUTPUT ReactionTurnOff_test_input_files
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/ReactionTurnOff/init_y .
  COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/tests/ReactionTurnOff/final_y .
)

add_executable(ReactionTurnOff ReactionTurnOff.cpp ReactionTurnOff_test_input_files)

add_test(ReactionTurnOff ReactionTurnOff)

target_link_libraries(ReactionTurnOff
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
