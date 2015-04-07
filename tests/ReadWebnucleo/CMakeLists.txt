add_executable(ReadWebnucleo ReadWebnucleo.cpp)

add_test(ReadWebnucleo ReadWebnucleo)

target_link_libraries(ReadWebnucleo
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
