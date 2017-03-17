/// \file Neutrino.hpp
/// \author lroberts
/// \since May 13, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_NEUTRINO_HPP_
#define SKYNET_REACTIONS_NEUTRINO_HPP_

#include <string>
#include <vector>

#include <H5Cpp.h>

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/NeutrinoEntry.hpp"
#include "Reactions/Reaction.hpp"

class Neutrino {
 public:
  explicit Neutrino(const std::string& pathToNeutrinoFile,
      const NuclideLibrary& nuclib);
   
  Neutrino(const std::vector<NeutrinoEntry>& entries, const std::string& source) 
      : mEntries(entries), mSource(source) {}
       
  static Neutrino FromREACLIBDecays(const std::string& REACLIBfile, 
      const NuclideLibrary& nuclib, NetworkOptions opt);
  
  std::vector<Reaction> GetValidReactions(
      const NuclideLibrary& nuclib) const;

  const std::vector<NeutrinoEntry>& Entries() const {
    return mEntries;
  }

  const std::vector<Reaction> Reactions(const NuclideLibrary& nuclib) const {
    std::vector<Reaction> reactions;
    for (auto entry : mEntries)
      reactions.push_back(entry.GetReaction(nuclib));

    return reactions;
  }

  std::string GetSource() const {
    return mSource;
  }

 private:
  std::vector<NeutrinoEntry> mEntries;
  std::string mSource;
};

#endif // SKYNET_REACTIONS_NEUTRINO_HPP_
