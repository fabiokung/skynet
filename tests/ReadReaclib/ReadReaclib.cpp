/// \file Reactlib.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <vector>

#include "BuildInfo.hpp"
#include "Reactions/REACLIBEntry.hpp"
#include "Reactions/REACLIB.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  if (REACLIB(SkyNetRoot + "/data/reaclib").Entries().size() != 82855)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
