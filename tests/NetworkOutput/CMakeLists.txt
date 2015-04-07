add_executable(NetworkOutput NetworkOutput.cpp)

add_test(NetworkOutput NetworkOutput)

target_link_libraries(NetworkOutput
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
