set(NuclearData_SOURCES
  Nuclide.cpp
  NuclideLibrary.cpp
  WebnucleoXmlReader.cpp
  WinvFileReader.cpp
)

add_SkyNet_library(NuclearData "${NuclearData_SOURCES}")

set(NuclearData_SWIG_DEPS
  Nuclide.i
  Nuclide.hpp
  NuclideLibrary.i
  NuclideLibrary.hpp
)

add_SWIG_dependencies("${NuclearData_SWIG_DEPS}")
