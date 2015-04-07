add_executable(NSE_screening NSE_screening.cpp)

add_test(NSE_screening NSE_screening)

target_link_libraries(NSE_screening
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
