/// \file REACLIB.hpp
/// \author jlippuner
/// \since Jul 4, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACLIB_HPP_
#define SKYNET_REACTIONS_REACLIB_HPP_

#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/REACLIBEntry.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name LeptonMode in Python instead of the enum values
// being global constants

#ifdef SWIG
%rename(LeptonMode) LeptonModeStruct;
#endif // SWIG

struct LeptonModeStruct {
  enum Value {
    TreatAllAsDecay,
    TreatAllAsCapture,
    TreatAllAsDecayExceptLabelEC
  };
};

typedef LeptonModeStruct::Value LeptonMode;

class REACLIB {
public:
  explicit REACLIB(const std::string& pathToREACLIBFile);

  static bool ReactionContainsUnkownNuclides(const REACLIBEntry& entry,
      const NuclideLibrary& nuclideLibrary);

  std::vector<Reaction> GetReactions(const ReactionType type,
      const NuclideLibrary& nuclideLibrary, const LeptonMode leptonMode) const;

  std::vector<REACLIBCoefficients_t> GetCoefficients(const ReactionType type,
      const NuclideLibrary& nuclideLibrary) const;

  std::vector<std::string> GetNuclideNames() const;

  const std::vector<REACLIBEntry>& Entries() const {
    return mEntries;
  }

private:
  REACLIBEntry ReadEntry(std::istream& is) const;

  void ReadFile(const std::string& pathToREACLIBFile);

  std::vector<REACLIBEntry> mEntries;
  std::unordered_set<std::string> mNuclideNames;
};

#endif // SKYNET_REACTIONS_REACLIB_HPP_
