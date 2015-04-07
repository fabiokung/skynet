set(FortranInterface_SOURCES
  FortranInterface.cpp
  SkyNet.f90
)

add_SkyNet_library(FortranInterface "${FortranInterface_SOURCES}")

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/skynet.mod
  DESTINATION "${INSTALL_INCLUDE_DIR}"
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
)
