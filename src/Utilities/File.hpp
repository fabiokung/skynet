/// \file File.hpp
/// \author jlippuner
/// \since Aug 28, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_FILE_HPP_
#define SKYNET_UTILITIES_FILE_HPP_

#include <string>

class File {
public:
  static bool Exists(const std::string path);

  static bool Writable(const std::string path);

  static std::string MakeAbsolutePath(const std::string path);
};

#endif // SKYNET_UTILITIES_FILE_HPP_
