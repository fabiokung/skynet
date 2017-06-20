add_executable(ReadFFN ReadFFN.cpp)

add_test(ReadFFN ReadFFN)

target_link_libraries(ReadFFN
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
