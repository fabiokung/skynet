add_executable(StopWatch StopWatch.cpp)

add_test(StopWatch StopWatch)

target_link_libraries(StopWatch
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
