/// \file WinvFileReader.hpp
/// \author jlippuner
/// \since Jul 28, 2014
///
/// \brief
///
///

#ifndef SKYNET_NUCLEARDATA_WINVFILEREADER_HPP_
#define SKYNET_NUCLEARDATA_WINVFILEREADER_HPP_

#include <string>
#include <unordered_map>

#include "NuclearData/Nuclide.hpp"

class WinvFileReader {
public:
  static std::unordered_map<std::string, Nuclide> GetAllNuclides(
      const std::string& pathToWinvFile);

private:
  // don't allow this class to be constructed
  WinvFileReader();
};

#endif // SKYNET_NUCLEARDATA_WINVFILEREADER_HPP_
