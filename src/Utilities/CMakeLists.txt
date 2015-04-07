add_subdirectory(Interpolators)

set(Utilities_SOURCES
  BisectionMethod.cpp
  CodeError.cpp
  File.cpp
  FloatingPointComparison.cpp
  FloatingPointExceptions.cpp
  FlushStdout.cpp
  FunctionIntegrator.cpp
  FunctionVsTimeWrapper.cpp
  #LineSearch.cpp
  OutputMuter.cpp
  Profiler.cpp
  StopWatch.cpp
)

add_SkyNet_library(Utilities "${Utilities_SOURCES}")

set(Utilities_SWIG_DEPS
  Constants.i
  Constants.hpp
  FunctionVsTime.i
  FunctionVsTime.hpp
)

add_SWIG_dependencies("${Utilities_SWIG_DEPS}")
