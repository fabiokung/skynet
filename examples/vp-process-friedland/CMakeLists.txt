add_executable(vp_process vp_process.cpp)

target_include_directories(vp_process PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(vp_process
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)

install(TARGETS vp_process
  RUNTIME DESTINATION "${INSTALL_EXAMPLES_DIR}"
)

install(DIRECTORY trajectories/
  DESTINATION "${INSTALL_EXAMPLES_DIR}/vp-process-friedland"
  OPTIONAL
)
