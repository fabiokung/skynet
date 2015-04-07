/// \file Reactlib.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <vector>

#include "BuildInfo.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  REACLIBReactionLibrary symmetricFission(
      SkyNetRoot + "/data/netsu_panov_symmetric_2neut",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Symmetric neutron induced fission with 2 free neutrons", nuclib,
      NetworkOptions(), false);
  REACLIBReactionLibrary neutronFission(
      SkyNetRoot + "/data/netsu_nfis_Roberts2010rates",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Neutron induced fission", nuclib, NetworkOptions(), false);
  REACLIBReactionLibrary spontaneousFission(
      SkyNetRoot + "/data/netsu_sfis_Roberts2010rates",
      ReactionType::Strong, false, LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Spontaneous fission", nuclib, NetworkOptions(), false);

  if (symmetricFission.NumAllReactions() != 3324) {
    printf("Wrong number of symmetric fission reactions: %lu instead of %i\n",
        symmetricFission.NumAllReactions(), 3324);
    return EXIT_FAILURE;
  }

  if (neutronFission.NumAllReactions() != 6309) {
    printf("Wrong number of neutron induced fission reactions: %lu instead of "
        "%i\n", neutronFission.NumAllReactions(), 6309);
    return EXIT_FAILURE;
  }

  if (spontaneousFission.NumAllReactions() != 10037) {
    printf("Wrong number of spontaneous fission reactions: %lu instead of %i\n",
        spontaneousFission.NumAllReactions(), 10037);
    return EXIT_FAILURE;
  }

  // spontaneous fissions must have a constant rate
  if (spontaneousFission.NumActiveNonConstReactions() != 0) {
    printf("spontaneous fission library has %lu non-constant reactions out of "
        "%lu\n", spontaneousFission.NumActiveNonConstReactions(),
        spontaneousFission.NumAllReactions());
    return EXIT_FAILURE;
  }

  // check the symmetric fissions make 2 neutrons and no protons
  int neutronId = nuclib.NuclideIdsVsNames().at("n");
  int protonId = nuclib.NuclideIdsVsNames().at("p");

  for (auto reaction : symmetricFission.Reactions().ActiveData()) {
    int neutronNum = 0;
    for (unsigned int i = 0; i < reaction.ProductNames().size(); ++i) {
      int id = nuclib.NuclideIdsVsNames().at(reaction.ProductNames()[i]);
      if (id == neutronId)
        neutronNum = reaction.NsOfProducts()[i];
      if (id == protonId) {
        printf("symmetric fission reaction contains proton product: %s\n",
            reaction.String().c_str());
        return EXIT_FAILURE;
      }
    }

    if (neutronNum != 2) {
      printf("symmetric fission reaction does not produce 2 neutrons: %s\n",
          reaction.String().c_str());
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
