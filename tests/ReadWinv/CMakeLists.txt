add_executable(ReadWinv ReadWinv.cpp)

add_test(ReadWinv ReadWinv)

target_link_libraries(ReadWinv
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
