/// \file FFN.cpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#include "Reactions/FFN.hpp"

#include "Utilities/File.hpp"

FFN::FFN(const std::string& pathToFile, const NuclideLibrary& nuclib) {
  // Open the the hdf5 file
  if (!File::Exists(pathToFile))
    throw std::invalid_argument(
        "The FFN file '" + pathToFile + "' does not exist.");

  mSource = File::MakeAbsolutePath(pathToFile);

  H5::H5File file(pathToFile.c_str(), H5F_ACC_RDONLY);

  // get the names of all the groups containing weak interactions
  std::vector<std::string> groupNames;
  groupNames.reserve(file.getNumObjs());
  for (unsigned int i = 0; i < file.getNumObjs(); ++i) {
    if (file.getObjTypeByIdx(i) == H5G_obj_t::H5G_GROUP) {
      char buf[1024];
      if (file.getObjnameByIdx(i, buf, 1024) == 0)
        throw std::invalid_argument(
          "Group without name in '" + pathToFile + "'");
      groupNames.push_back(std::string(buf));
    }
  }

  // load the reactions
  mEntries = std::vector<FFNEntry>(2 * groupNames.size());
  mReactions = std::vector<Reaction>(2 * groupNames.size());

  for (unsigned int i = 0; i < groupNames.size(); ++i) {
    auto entry = FFNEntry(file, groupNames[i]);
    mEntries[2 * i] = entry;
    mEntries[2 * i + 1] = entry;

    auto reacs = entry.GetReactions(nuclib);
    mReactions[2 * i] = reacs[0];
    mReactions[2 * i + 1] = reacs[1];
  }

  if (mEntries.size() > 0) {
    mT9Grid = mEntries[0].T9Grid();
    mLogYeRhoGrid = mEntries[0].LogYeRhoGrid();
  }

  // perform some simple error checking
  for (auto entry : mEntries) {
    if (entry.T9Grid() != mT9Grid)
      throw std::runtime_error(entry.Label() + " has different T9Grid");
    if (entry.LogYeRhoGrid() != mLogYeRhoGrid)
      throw std::runtime_error(entry.Label() + " has different LogYeRhoGrid");
  }

  if (!std::is_sorted(mT9Grid.begin(), mT9Grid.end()))
    throw std::runtime_error("Non-monotonically increasing T9 grid in FFN "
        "reactions");

  if (!std::is_sorted(mLogYeRhoGrid.begin(), mLogYeRhoGrid.end()))
    throw std::runtime_error("Non-monotonically increasing logYeRho grid in "
        "FFN reactions");
}

std::vector<Reaction> FFN::GetValidReactions(
    const NuclideLibrary& nuclib) const {
  std::vector<Reaction> reactions;
  for (auto reaction : mReactions) {
    if (reaction.ContainsOnlyKnownNuclides(nuclib))
      reactions.push_back(reaction);
  }
  return reactions;
}

// get vector of one point in temperature density space table for all reactions
std::vector<double> FFN::GetTableEntryAll(const int YeRhoIdx,
    const int T9Idx, const NuclideLibrary& nuclib) const {
  std::vector<double> data;
  int idx = GetTableIdx(YeRhoIdx, T9Idx);
  for (unsigned int i = 0; i < mEntries.size(); ++i) {
    if (mReactions[i].ContainsOnlyKnownNuclides(nuclib)) {
      if (mReactions[i].IsDecay())
        data.push_back(mEntries[i].Beta()[idx]);
      else
        data.push_back(mEntries[i].Eps()[idx]);
    }
  }
  return data;
}

// Write rate table to screen for testing purposes
//void FFN::writeReaction(int idx) const
//    {
//  if ((unsigned int)idx < mReactions.size()) {
//    const FFNReaction *reac = &mEntries[idx];
//    int wd = 8;
//    cout << std::string(6 * wd, ' ') << reac->nameParent << " -> "
//        << reac->nameDaughter << endl;
//    cout << std::setw(wd + 3) << " ";
//    for (auto T9 : reac->T9Grid)
//      cout << std::setw(wd) << T9 << " ";
//    cout << endl;
//    cout << std::setw(wd + 3) << " ";
//    for (unsigned int i = 0; i < reac->T9Grid.size(); i++)
//      cout << std::string(wd + 1, '-');
//    cout << endl;
//    for (unsigned int i = 0; i < reac->logYeRhoGrid.size(); i++) {
//      cout << std::setw(wd) << reac->logYeRhoGrid[i] << " | ";
//      for (unsigned int j = 0; j < reac->T9Grid.size(); j++)
//        cout << std::setw(wd) << reac->beta[reac->getTabIdx(i, j)] << " ";
//      cout << endl;
//    }
//  }
//  else {
//    cout << " That reaction isn't available." << endl;
//  }
//}
