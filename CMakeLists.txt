cmake_minimum_required(VERSION 3.5)

enable_testing()

# set Release build type as default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Configuration type (one of Debug, RelWithDebInfo, Release, MinSizeRel)"
    FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules/")

project(SkyNet CXX C Fortran)

# enable warnings
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# C++14 required by Armadillo >= 10 (original code was C++11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

# generate debugging symbols for release and debug
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")

# explicitly set DEBUG flag in Debug mode
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

# to make shared libraries we need position independent code
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

# linking the shared Intel libraries suppresses warnings about feupdateenv
# not being implemented, and to be sure that we're not missing any math functions
# we also link against libm
if(("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel"))
  set(SKYNET_EXTERNAL_LIBS "imf;m;${SKYNET_EXTERNAL_LIBS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -shared-intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -shared-intel")
endif ()

# if SKYNET_ROOT does not exist or if the CMAKE_INSTALL_PREFIX is different from
# the cached one, set SKYNET_ROOT and cache it
if((NOT SKYNET_ROOT) OR (NOT (CMAKE_INSTALL_PREFIX STREQUAL CMAKE_INSTALL_PREFIX_CACHED)))
  if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(SKYNET_ROOT "${CMAKE_SOURCE_DIR}" CACHE INTERNAL "SkyNetRoot")
  else()
    set(SKYNET_ROOT "${CMAKE_INSTALL_PREFIX}" CACHE INTERNAL "SkyNetRoot")
  endif()
  set(CMAKE_INSTALL_PREFIX_CACHED "${CMAKE_INSTALL_PREFIX}" CACHE INTERNAL "Cache of CMAKE_INSTALL_PREFIX")
endif()

# offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_EXAMPLES_DIR examples CACHE PATH "Installation directory for examples")
set(INSTALL_DATA_DIR data CACHE PATH "Installation directory for data")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/cmake/SkyNet)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE EXAMPLES DATA CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

include("CMakeLists_Requirements.txt")
include("CMakeLists_MatrixSolver.txt")

# add the source directory to the include list
include_directories(${SkyNet_SOURCE_DIR}/src)

macro(add_SkyNet_library library sources)
  add_library(${library} OBJECT ${sources})
  set(SKYNET_LIBS "${SKYNET_LIBS};$<TARGET_OBJECTS:${library}>" PARENT_SCOPE)
  add_dependencies(${library} src)
endmacro()

macro(add_SWIG_dependencies)
  foreach(dep ${ARGN})
    list(APPEND SWIG_MODULE_SkyNet_EXTRA_DEPS "${CMAKE_CURRENT_SOURCE_DIR}/${dep}")
  endforeach(dep)
  set(SWIG_MODULE_SkyNet_EXTRA_DEPS ${SWIG_MODULE_SkyNet_EXTRA_DEPS} PARENT_SCOPE)
endmacro()

add_subdirectory(src)

# create static library that will be used to link tests and SWIG library
add_library(SkyNet_static
  ${SKYNET_LIBS}
)

# create shared library that can be used to link other C++ applications
add_library(SkyNet SHARED
  ${SKYNET_LIBS}
)

# create SWIG Python bindings
if(USE_SWIG)
  SET_SOURCE_FILES_PROPERTIES("src/SkyNet.i" PROPERTIES CPLUSPLUS ON)
  if(USE_MOVIE)
    set_property(SOURCE "src/SkyNet.i" PROPERTY SWIG_FLAGS "-DSKYNET_USE_MOVIE")
  endif()

  if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 3.8)
    SWIG_ADD_MODULE(SkyNet python "src/SkyNet.i")
  else()
    SWIG_ADD_LIBRARY(SkyNet LANGUAGE python SOURCES "src/SkyNet.i")
  endif()

  SWIG_LINK_LIBRARIES(SkyNet
    SkyNet_static
    ${SKYNET_EXTERNAL_LIBS}
    ${Python3_LIBRARIES}
  )

  install(FILES ${CMAKE_BINARY_DIR}/_SkyNet.so
    DESTINATION "${INSTALL_LIB_DIR}"
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )

  install(FILES ${CMAKE_BINARY_DIR}/SkyNet.py
    DESTINATION "${INSTALL_LIB_DIR}"
  )
endif()

add_subdirectory(tests)
add_subdirectory(examples)

# define install target
install(TARGETS SkyNet
  EXPORT SkyNetTargets
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

install(TARGETS SkyNet_static
  EXPORT SkyNetTargets
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
)

install(DIRECTORY data/
  DESTINATION "${INSTALL_DATA_DIR}"
  USE_SOURCE_PERMISSIONS
)

# make CMake config file

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
set(SKYNET_INCLUDE_DIRS "\${SKYNET_CMAKE_DIR}/${REL_INCLUDE_DIR}")

# get external include directories
get_target_property(EXTERNAL_INCLUDE_DIRS SkyNet INCLUDE_DIRECTORIES)
list(APPEND SKYNET_INCLUDE_DIRS "${EXTERNAL_INCLUDE_DIRS}")

configure_file(SkyNetConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SkyNetConfig.cmake" @ONLY
)

# install the config file
install(FILES "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/SkyNetConfig.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}"
)

# install the export set for use with the install-tree
install(EXPORT SkyNetTargets
  DESTINATION "${INSTALL_CMAKE_DIR}"
)

# the compiler flags
message(STATUS "")
message(STATUS "")
message(STATUS "Configuration options:")
message(STATUS "")
message(STATUS "              CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "         CMAKE_CXX_LINK_FLAGS: " ${CMAKE_CXX_LINK_FLAGS})
message(STATUS "        CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "      CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "          CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})
message(STATUS "     CMAKE_Fortran_LINK_FLAGS: " ${CMAKE_Fortran_LINK_FLAGS})
message(STATUS "    CMAKE_Fortran_FLAGS_DEBUG: " ${CMAKE_Fortran_FLAGS_DEBUG})
message(STATUS "  CMAKE_Fortran_FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})

message(STATUS "")
message(STATUS "")

if(SKYNET_USE_DENSE_SOLVER)
  message(STATUS "                Matrix solver: " ${SKYNET_MATRIX_SOLVER} " (dense)")
else()
  message(STATUS "                Matrix solver: " ${SKYNET_MATRIX_SOLVER} " (sparse)")
endif()

if(USE_SWIG)
  message(STATUS "                         SWIG: enabled")
else()
  message(STATUS "                         SWIG: disabled")
endif()

if(USE_MOVIE)
  message(STATUS "                making movies: enabled")
else()
  message(STATUS "                making movies: disabled")
endif()

message(STATUS "")
message(STATUS "")
