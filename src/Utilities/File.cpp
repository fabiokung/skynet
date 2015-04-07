/// \file File.cpp
/// \author jlippuner
/// \since Aug 28, 2014
///
/// \brief
///
///

#include "Utilities/File.hpp"

#include <cstdio>

#include <boost/filesystem.hpp>

bool File::Exists(const std::string path) {
  FILE * f = fopen(path.c_str(), "r");
  if (f != nullptr) {
    fclose(f);
    return true;
  } else {
    return false;
  }
}

bool File::Writable(const std::string path) {
  FILE * f = fopen(path.c_str(), "a");
  if (f != nullptr) {
    fclose(f);
    return true;
  } else {
    return false;
  }
}

std::string File::MakeAbsolutePath(const std::string path) {
  // if source is a file, get its aboslute path
  boost::filesystem::path file(path);

  if (boost::filesystem::exists(file))
    return boost::filesystem::absolute(file).string();
  else
    return path;
}
