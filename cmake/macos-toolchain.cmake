# macOS toolchain: Apple Clang with Command Line Tools
#
# On macOS the CLT ships C++ headers in the SDK only (not the standalone
# include path), so we point the compiler at the SDK explicitly.
#
# The shared library build also needs -undefined dynamic_lookup because
# CMakeLists.txt does not call target_link_libraries on the shared target.

execute_process(
    COMMAND xcrun --show-sdk-path
    OUTPUT_VARIABLE _macos_sdk_path
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

if(_macos_sdk_path)
    set(CMAKE_OSX_SYSROOT "${_macos_sdk_path}" CACHE PATH "macOS SDK root")
    set(CMAKE_CXX_FLAGS_INIT "-isystem ${_macos_sdk_path}/usr/include/c++/v1")
    set(CMAKE_SHARED_LINKER_FLAGS_INIT "-undefined dynamic_lookup")
endif()
