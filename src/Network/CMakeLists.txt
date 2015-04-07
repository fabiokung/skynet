set(Network_SOURCES
  NetworkOptions.cpp
  NetworkOutput.cpp
  NetworkOutputDataset.cpp
  NSE.cpp
  ReactionNetwork.cpp
  TemperatureDensityHistory.cpp
)

add_SkyNet_library(Network "${Network_SOURCES}")

set(Network_SWIG_DEPS
  NetworkOptions.i
  NetworkOptions.hpp
  NetworkOutput.i
  NetworkOutput.hpp
  NSE.i
  NSE.hpp
  NSEOptions.i
  NSEOptions.hpp
  ReactionNetwork.i
  ReactionNetwork.hpp
  TemperatureDensityHistory.i
  TemperatureDensityHistory.hpp
)

add_SWIG_dependencies("${Network_SWIG_DEPS}")

