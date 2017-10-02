add_executable(InternalEnergyConservation InternalEnergyConservation.cpp)

add_test(InternalEnergyConservation InternalEnergyConservation)

target_link_libraries(InternalEnergyConservation
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
