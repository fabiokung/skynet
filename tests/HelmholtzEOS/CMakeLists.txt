add_executable(HelmholtzEOSTest HelmholtzEOS.cpp)

add_test(HelmholtzEOS HelmholtzEOSTest)

target_link_libraries(HelmholtzEOSTest
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
