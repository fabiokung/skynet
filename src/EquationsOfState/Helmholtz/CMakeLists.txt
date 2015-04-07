set(HelmholtzEOS_SOURCES
  helmholtz.f90
  helmholtz_addition.f90
)

add_SkyNet_library(HelmholtzEOS "${HelmholtzEOS_SOURCES}")

install(FILES helm_table.dat
  DESTINATION "${INSTALL_DATA_DIR}"
)
