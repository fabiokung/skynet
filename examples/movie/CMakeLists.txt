add_executable(movie movie.cpp)

target_link_libraries(movie
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)

install(TARGETS movie
  RUNTIME DESTINATION examples
)
