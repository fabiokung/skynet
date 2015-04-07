set(MovieMaker_SOURCES
  Scale.cpp
  ColorScale.cpp
  PiecewiseLinearColorScale.cpp
  CairoImage.cpp
  NucleiChart.cpp
  MovieData.cpp
)

add_SkyNet_library(MovieMaker "${MovieMaker_SOURCES}")

set(MovieMaker_SWIG_DEPS
  ColorScale.i
  ColorScale.hpp
  ColorScales.i
  ColorScales.hpp
  MovieData.i
  MovieData.hpp
  NucleiChart.i
  NucleiChart.hpp
  NucleiChartOptions.i
  NucleiChartOptions.hpp
  PiecewiseLinearColorScale.i
  PiecewiseLinearColorScale.hpp
  Scale.i
  Scale.hpp
)

add_SWIG_dependencies("${MovieMaker_SWIG_DEPS}")
