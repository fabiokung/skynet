add_executable(InverseRates InverseRates.cpp)

add_test(InverseRates InverseRates)

target_link_libraries(InverseRates
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
