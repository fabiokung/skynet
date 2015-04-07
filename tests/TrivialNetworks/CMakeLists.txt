add_executable(TrivialNetworks TrivialNetworks.cpp )

add_test(TrivialNetworks TrivialNetworks)

target_link_libraries(TrivialNetworks
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
