add_executable(NeutrinoNetwork NeutrinoNetwork.cpp)

add_test(NeutrinoNetwork NeutrinoNetwork)

target_link_libraries(NeutrinoNetwork
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
