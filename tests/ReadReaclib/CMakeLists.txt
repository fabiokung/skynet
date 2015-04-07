add_executable(ReadReaclib ReadReaclib.cpp)

add_test(ReadReaclib ReadReaclib)

target_link_libraries(ReadReaclib
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
