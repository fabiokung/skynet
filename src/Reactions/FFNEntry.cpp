/// \file FFNEntry.cpp
/// \author jlippuner
/// \since Jan 31, 2015
///
/// \brief
///
///

#include "Reactions/FFNEntry.hpp"

#include <algorithm>
#include <numeric>
#include <set>

#include "NuclearData/Nuclide.hpp"
#include "Utilities/H5Helper.hpp"

namespace { // unnamed so this function can only be used in this file

template<typename T>
std::vector<T> ReadData(const H5::Group& group, const std::string& name) {
  auto datSpace = group.openDataSet(name.c_str()).getSpace();

  hsize_t ndim;
  ndim = datSpace.getSimpleExtentNdims();

  hsize_t dims[ndim];
  ndim = datSpace.getSimpleExtentDims(dims, NULL);

  int ntot = 1;
  for (unsigned int i = 0; i < ndim; ++i)
    ntot *= dims[i];

  std::vector<T> dat(ntot);
  group.openDataSet(name.c_str()).read(dat.data(), GetH5DataType<T>());

  return dat;
}

}

FFNEntry::FFNEntry(const H5::H5File& h5File, const std::string& groupName) {
  try {
    H5::Group group = h5File.openGroup(groupName.c_str());

    // read in general reaction data
    mParentZs = ReadData<int>(group, "Parent_Z");
    mParentAs = ReadData<int>(group, "Parent_A");
    for (unsigned int i = 0; i < mParentZs.size(); i++)
      mParentNames.push_back(Nuclide::GetName(mParentZs[i], mParentAs[i]));

    mDaughterZs = ReadData<int>(group, "Daughter_Z");
    mDaughterAs = ReadData<int>(group, "Daughter_A");
    for (unsigned int i = 0; i < mDaughterZs.size(); i++)
      mDaughterNames.push_back(
          Nuclide::GetName(mDaughterZs[i], mDaughterAs[i]));

    // Create a human readable identifier for the reaction
    mLabel = "";
    for (auto name : mParentNames)
      mLabel += name + " ";

    std::string daughters = "";
    for (auto name : mDaughterNames)
      daughters += " " + name;

    if (std::accumulate(mParentZs.begin(), mParentZs.end(), 0)
        > std::accumulate(mDaughterZs.begin(), mDaughterZs.end(), 0))
      mLabel += "+ {  ,e-} ->" + daughters + " + {e+,  } + nue";
    else
      mLabel += "+ {  ,e+} ->" + daughters + " + {e-,  } + nueb";

    // read in temperature and density grids
    mT9Grid = ReadData<double>(group, "T9grid");
    if (mT9Grid.size() <= 0)
      throw std::invalid_argument("T9grid is empty");

    mLogYeRhoGrid = ReadData<double>(group, "logYeRhogrid");
    if (mLogYeRhoGrid.size() <= 0)
      throw std::invalid_argument("logYeRhogrid is empty");

    // get the rate array size
    mBeta = ReadData<double>(group, "beta");
    if (mBeta.size() != mT9Grid.size() * mLogYeRhoGrid.size())
      throw std::invalid_argument("beta has wrong size");

    mEps = ReadData<double>(group, "eps");
    if (mEps.size() != mT9Grid.size() * mLogYeRhoGrid.size())
      throw std::invalid_argument("eps has wrong size");

    mLnu = ReadData<double>(group, "lnu");
    if (mLnu.size() != mT9Grid.size() * mLogYeRhoGrid.size())
      throw std::invalid_argument("lnu has wrong size");
  } catch (std::exception& ex) {
    throw std::runtime_error("Error while reading FFN group '" + groupName
        + "', inner exception:\n" + ex.what());
  }
}

std::vector<Reaction> FFNEntry::GetReactions(
    const NuclideLibrary& nuclib) const {
  // Find unique reactants and products
  std::set<std::string> uniqueParents, uniqueDaughters;
  for (auto name : mParentNames)
    uniqueParents.insert(name);
  for (auto name : mDaughterNames)
    uniqueDaughters.insert(name);

  // Count up the number of each unique reactant and product
  std::vector<int> NReacs, NProds;
  for (auto name : uniqueParents)
    NReacs.push_back(
        std::count(mParentNames.begin(), mParentNames.end(), name));
  for (auto name : uniqueDaughters)
    NProds.push_back(
        std::count(mDaughterNames.begin(), mDaughterNames.end(), name));

  Reaction decay(
      std::vector<std::string>(uniqueParents.begin(), uniqueParents.end()),
      std::vector<std::string>(uniqueDaughters.begin(), uniqueDaughters.end()),
      NReacs, NProds, true, false, true, mLabel, nuclib);

  Reaction capture(
      std::vector<std::string>(uniqueParents.begin(), uniqueParents.end()),
      std::vector<std::string>(uniqueDaughters.begin(), uniqueDaughters.end()),
      NReacs, NProds, true, false, false, mLabel, nuclib);

  return {decay, capture};
}
