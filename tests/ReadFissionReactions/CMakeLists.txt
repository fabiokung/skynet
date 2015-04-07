add_executable(ReadFissionReactions ReadFissionReactions.cpp)

add_test(ReadFissionReactions ReadFissionReactions)

target_link_libraries(ReadFissionReactions
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
