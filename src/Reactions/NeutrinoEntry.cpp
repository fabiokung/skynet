/// \file NeutrinoEntry.cpp
/// \author lroberts
/// \since May 13, 2015
///
/// \brief
///
///

#include "Reactions/NeutrinoEntry.hpp"

#include <algorithm>
#include <numeric>
#include <set>

#include "NuclearData/Nuclide.hpp"
#include "Utilities/H5Helper.hpp"

NeutrinoEntry::NeutrinoEntry(const std::vector<int> ParentAs,
    const std::vector<int> ParentZs, const std::vector<int> DaughterAs,
    const std::vector<int> DaughterZs, const double Q,
    const double matrixElement, const double Wm, const bool isNue) :
    mQ(Q), mMatrixElement(matrixElement), mWm(Wm), mIsNue(isNue) {
  // Build the names
  for (unsigned int i = 0; i<ParentAs.size(); ++i)
     mParentNames.push_back(Nuclide::GetName(ParentZs[i], ParentAs[i]));
  for (unsigned int i = 0; i<DaughterAs.size(); ++i)
     mDaughterNames.push_back(Nuclide::GetName(DaughterZs[i], DaughterAs[i]));
}

Reaction NeutrinoEntry::GetReaction(const NuclideLibrary& nuclib) const {
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

  return Reaction(
      std::vector<std::string>(uniqueParents.begin(), uniqueParents.end()),
      std::vector<std::string>(uniqueDaughters.begin(), uniqueDaughters.end()),
      NReacs, NProds, true, false, false, "Neutrino Reaction", nuclib);
}
