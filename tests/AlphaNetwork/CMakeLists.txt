add_executable(AlphaNetwork AlphaNetwork.cpp)

add_test(AlphaNetwork AlphaNetwork)

target_link_libraries(AlphaNetwork
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
