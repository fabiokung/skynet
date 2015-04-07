/// \file ReadWebnucleo.cpp
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
#include "NuclearData/NuclideLibrary.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  NuclideLibrary nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  if (nuclib.NumNuclides() != 8230)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
