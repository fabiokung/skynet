add_executable(SmallNetworks SmallNetworks.cpp )

add_test(SmallNetworks SmallNetworks)

target_link_libraries(SmallNetworks
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
