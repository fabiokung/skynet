/// \file FFNEntry.hpp
/// \author jlippuner
/// \since Jan 31, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_FFNENTRY_HPP_
#define SKYNET_REACTIONS_FFNENTRY_HPP_

#include <string>
#include <vector>

#include <H5Cpp.h>

#include "Reactions/Reaction.hpp"

class FFNEntry {
public:
  FFNEntry() {}

  FFNEntry(const H5::H5File& h5File, const std::string& groupName);

  std::vector<Reaction> GetReactions(const NuclideLibrary& nuclib) const;

  const std::vector<std::string>& ParentNames() const {
    return mParentNames;
  }

  const std::vector<std::string>& NameDaughter() const {
    return mDaughterNames;
  }

  std::vector<int> ParentZs() const {
    return mParentZs;
  }

  std::vector<int> ParentAs() const {
    return mParentAs;
  }

  std::vector<int> DaughterZs() const {
    return mDaughterZs;
  }

  std::vector<int> DaughterAs() const {
    return mDaughterAs;
  }

  const std::string& Label() const {
    return mLabel;
  }

  const std::vector<double>& Beta() const {
    return mBeta;
  }

  const std::vector<double>& Eps() const {
    return mEps;
  }

  const std::vector<double>& Lnu() const {
    return mLnu;
  }

  const std::vector<double>& T9Grid() const {
    return mT9Grid;
  }

  const std::vector<double>& LogYeRhoGrid() const {
    return mLogYeRhoGrid;
  }

private:
  std::vector<std::string> mParentNames, mDaughterNames;

  std::vector<int> mParentZs, mParentAs;
  std::vector<int> mDaughterZs, mDaughterAs;

  std::string mLabel;

  std::vector<double> mBeta;
  std::vector<double> mEps;
  std::vector<double> mLnu;

  std::vector<double> mT9Grid;
  std::vector<double> mLogYeRhoGrid;
};

#endif // SKYNET_REACTIONS_FFNENTRY_HPP_
