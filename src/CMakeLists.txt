add_subdirectory(DensityProfiles)
add_subdirectory(EquationsOfState)
add_subdirectory(FortranInterface)
add_subdirectory(MatrixSolver)
if(USE_MOVIE)
  add_subdirectory(MovieMaker)
endif()
add_subdirectory(Network)
add_subdirectory(NuclearData)
add_subdirectory(Reactions)
add_subdirectory(Utilities)
add_subdirectory(tinyxml2)

# get the Git root directory
execute_process(
  COMMAND git rev-parse --show-toplevel
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  RESULT_VARIABLE CMAKE_GIT_RESULT
  OUTPUT_VARIABLE CMAKE_GIT_ROOT
  ERROR_QUIET
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(CMAKE_GIT_RESULT EQUAL 0)
  set(CMAKE_USE_GIT "true")
else()
  set(CMAKE_USE_GIT "false")
  set(CMAKE_GIT_ROOT "")
endif()

# this creates MakeBuildInfo.cmake with the correct paths
configure_file(
  ${CMAKE_SOURCE_DIR}/CMakeModules/MakeBuildInfo.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MakeBuildInfo.cmake @ONLY
)

# this will run MakeBuildInfo.cmake at every make to update the time and Git
# hash
add_custom_target(MakeBuildInfo ALL
  COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/MakeBuildInfo.cmake
)

# this runs MakeBuildInfo.cmake in the configuration phase, so that the
# input file BuildInfo.cpp is created
include(${CMAKE_CURRENT_BINARY_DIR}/MakeBuildInfo.cmake)

set(src_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp
)

add_SkyNet_library(src "${src_SOURCES}")
add_dependencies(src MakeBuildInfo)

set(src_SWIG_DEPS
  BuildInfo.i
  BuildInfo.hpp
  SkyNet.i
  std_array.i
)

add_SWIG_dependencies("${src_SWIG_DEPS}")

install(DIRECTORY ./
  DESTINATION "${INSTALL_INCLUDE_DIR}"
  FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
)
