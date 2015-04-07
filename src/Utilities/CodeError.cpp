#include "Utilities/CodeError.hpp"

#include <iostream>

void AssertFail(const char* expression, const char* file, const int line,
    std::string msg) {
  std::cerr << std::endl
      << "#######################  ASSERT  FAILED  #####################"
      << std::endl
      << "#### " << expression << "  violated"
      << std::endl
      << "#### " << file << " (line " << line << ")"
      << std::endl
      << msg
      << std::endl
      << "###############################################################"
      << std::endl;
  throw CodeError();
}
