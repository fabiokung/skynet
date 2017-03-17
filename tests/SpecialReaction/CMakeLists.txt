add_executable(SpecialReaction SpecialReaction.cpp)

add_test(SpecialReaction SpecialReaction)

target_link_libraries(SpecialReaction
  SkyNet_static
  ${SKYNET_EXTERNAL_LIBS}
)
