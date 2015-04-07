# - Try to find Pardiso library (http://www.pardiso-project.org/)
# Once done, this will define
#
#  Pardiso_FOUND - system has Cairomm
#  Pardiso_LIBRARIES - link these to use Pardiso

include(LibFindMacros)

# find the library
find_library(Pardiso_LIBRARY
  NAMES pardiso pardiso500-GNU481-X86-64 pardiso500-INTEL1301-X86-64
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries that this lib depends on.
set(Pardiso_PROCESS_LIBS Pardiso_LIBRARY Pardiso_LIBRARIES)
libfind_process(Pardiso)
