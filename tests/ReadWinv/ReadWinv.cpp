/// \file ReadWinv.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <string>
#include <unordered_map>
#include <vector>

#include "BuildInfo.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/WinvFileReader.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  std::unordered_map<std::string, Nuclide> nuclides =
      WinvFileReader::GetAllNuclides(SkyNetRoot + "/data/winvn_v2.0.dat");

  if (nuclides.size() != 7852)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
