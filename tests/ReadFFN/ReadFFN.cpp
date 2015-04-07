/// \file ReadFFN.cpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#include <vector>

#include "BuildInfo.hpp"
#include "Reactions/FFN.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  FFN lib("mesa_weak_rates.h5", nuclib);
  std::vector<Reaction> reacs = lib.Reactions();

  for (auto reac : reacs)
    printf("%s\n", reac.Label().c_str());
  //lib.WriteReaction(1);

  return EXIT_SUCCESS;
}
