set(Utilities_Interpolators_SOURCES
  CubicHermiteInterpolator.cpp
  PiecewiseLinearFunction.cpp
)

add_SkyNet_library(Utilities_Interpolators "${Utilities_Interpolators_SOURCES}")

set(Utilities_Interpolators_SWIG_DEPS
  CubicHermiteInterpolator.i
  CubicHermiteInterpolator.hpp
  PiecewiseLinearFunction.i
  PiecewiseLinearFunction.hpp
  InterpolatorBase.i
  InterpolatorBase.hpp
  InterpolatorBase_defns.hpp
)

add_SWIG_dependencies("${Utilities_Interpolators_SWIG_DEPS}")
