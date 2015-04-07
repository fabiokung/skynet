set(DensityProfiles_SOURCES
  ConstantFunction.cpp
  ExpTMinus3.cpp
  PowerLawContinuation.cpp
)

add_SkyNet_library(DensityProfiles "${DensityProfiles_SOURCES}")

set(DensityProfiles_SWIG_DEPS
  ConstantFunction.i
  ConstantFunction.hpp
  ExpTMinus3.i
  ExpTMinus3.hpp
  PowerLawContinuation.i
  PowerLawContinuation.hpp
)

add_SWIG_dependencies("${DensityProfiles_SWIG_DEPS}")
