add_executable(r-process r-process.cpp)

target_link_libraries(r-process
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)

install(TARGETS r-process
  RUNTIME DESTINATION "${INSTALL_EXAMPLES_DIR}"
)
