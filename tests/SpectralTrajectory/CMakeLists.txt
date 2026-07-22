add_executable(SpectralTrajectory SpectralTrajectory.cpp)

target_include_directories(SpectralTrajectory PRIVATE
  ${CMAKE_SOURCE_DIR}/examples/vp-process-friedland)

add_test(SpectralTrajectory SpectralTrajectory)

target_link_libraries(SpectralTrajectory
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
