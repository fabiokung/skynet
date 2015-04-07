/// \file FFN.hpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_FFN_HPP_
#define SKYNET_REACTIONS_FFN_HPP_

#include <string>
#include <vector>

#include <H5Cpp.h>

#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/FFNEntry.hpp"
#include "Reactions/Reaction.hpp"

class FFN {
public:
  explicit FFN(const std::string& pathToFFNFile, const NuclideLibrary& nuclib);

  std::vector<Reaction> GetValidReactions(const NuclideLibrary& nuclib) const;

  std::vector<double> GetTableEntryAll(const int YeRhoIdx, const int T9Idx,
      const NuclideLibrary& nuclib) const;

//  void WriteReaction(const int idx) const;

  const std::vector<double>& T9Grid() const {
    return mT9Grid;
  }

  const std::vector<double>& LogYeRhoGrid() const {
    return mLogYeRhoGrid;
  }

  const std::vector<FFNEntry>& Entries() const {
    return mEntries;
  }

  const std::vector<Reaction>& Reactions() const {
    return mReactions;
  }

  std::string getSource() const {
    return mSource;
  }

private:
  inline int GetTableIdx(const int YeRhoIdx, const int T9Idx) const {
    return T9Idx + mT9Grid.size() * YeRhoIdx;
  }

  std::vector<double> mT9Grid;
  std::vector<double> mLogYeRhoGrid;

  std::vector<FFNEntry> mEntries;
  std::vector<Reaction> mReactions;

  std::string mSource;
};

#endif // SKYNET_REACTIONS_FFN_HPP_
