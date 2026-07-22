add_executable(TabulatedNeutrinoSpectrum TabulatedNeutrinoSpectrum.cpp)

add_test(TabulatedNeutrinoSpectrum TabulatedNeutrinoSpectrum)

target_link_libraries(TabulatedNeutrinoSpectrum
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
