add_executable(FunctionIntegrator_test FunctionIntegrator_test.cpp)

add_test(FunctionIntegrator FunctionIntegrator_test)

target_link_libraries(FunctionIntegrator_test
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
