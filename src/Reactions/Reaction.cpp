/// \file Reaction.cpp
/// \author jlippuner
/// \since Jun 30, 2014
///
/// \brief
///
///

#include "Reactions/Reaction.hpp"

#include <algorithm>
#include <sstream>

namespace { // unnamed so this function can only be used in this file

void SortNamesAndNumbers(std::vector<std::string> * const pNames,
    std::vector<int> * const pNumbers) {
  std::vector<std::pair<std::string, int>> vec(pNames->size());

  for (unsigned int i = 0; i < pNames->size(); ++i)
    vec[i] = { (*pNames)[i], (*pNumbers)[i] };

  std::sort(vec.begin(), vec.end(),
      [] (const std::pair<std::string, int> first,
          const std::pair<std::string, int> second) {
        return first.first < second.first;
      });

  for (unsigned int i = 0; i < pNames->size(); ++i) {
    (*pNames)[i] = vec[i].first;
    (*pNumbers)[i] = vec[i].second;
  }
}

}

std::string Leptons::ToString(const Lepton lepton) {
  switch (lepton) {
    case         Lepton::E: return "e-";
    case     Lepton::AntiE: return "e+";
    case        Lepton::Mu: return "mu-";
    case    Lepton::AntiMu: return "mu+";
    case       Lepton::Tau: return "tau-";
    case   Lepton::AntiTau: return "tau+";
    case       Lepton::NuE: return "nue";
    case   Lepton::AntiNuE: return "nueb";
    case       Lepton::NuX: return "nux";
    case      Lepton::NuMu: return "numu";
    case  Lepton::AntiNuMu: return "numub";
    case     Lepton::NuTau: return "nutau";
    case Lepton::AntiNuTau: return "nutaub";
    default:                return "<UNKNOWN LEPTON>";
  }
}

// (e, Le, Lmu, Ltau) where e is electric charge, Le is electronic lepton
// number, Lmu is muonic lepton number, and Ltau is tauonic lepton number
std::valarray<int> Leptons::Charge(const Lepton lepton) {
  switch (lepton) {
    case         Lepton::E: return { -1,  1,  0,  0 };
    case     Lepton::AntiE: return {  1, -1,  0,  0 };
    case        Lepton::Mu: return { -1,  0,  1,  0 };
    case    Lepton::AntiMu: return {  1,  0, -1,  0 };
    case       Lepton::Tau: return { -1,  0,  0,  1 };
    case   Lepton::AntiTau: return {  1,  0,  0, -1 };
    case       Lepton::NuE: return {  0,  1,  0,  0 };
    case   Lepton::AntiNuE: return {  0, -1,  0,  0 };
    case       Lepton::NuX: return {  0,  0,  0,  0 };
    case      Lepton::NuMu: return {  0,  0,  1,  0 };
    case  Lepton::AntiNuMu: return {  0,  0, -1,  0 };
    case     Lepton::NuTau: return {  0,  0,  0,  1 };
    case Lepton::AntiNuTau: return {  0,  0,  0, -1 };
    default:                return {  0,  0,  0,  0 };
  }
}

Reaction::Reaction(const std::vector<std::string>& reactantNames,
    const std::vector<std::string>& productNames,
    const std::vector<int>& NsOfReactants,
    const std::vector<int>& NsOfProducts, const bool isWeak,
    const bool isInverse, const bool isDecay, const std::string& label,
    const NuclideLibrary& nuclib) :
    mReactantNames(reactantNames),
    mProductNames(productNames),
    mNsOfReactants(NsOfReactants),
    mNsOfProducts(NsOfProducts),
    mIsWeak(isWeak),
    mIsInverse(isInverse),
    mLabel(label) {
  Init();

  // add leptons (only electron type)
  int chargeDiff = 0;
  for (unsigned int i = 0; i < mReactantNames.size(); ++i) {
    if (nuclib.NuclideIdsVsNames().count(mReactantNames[i]) == 0)
      // ignore this reaction, it will be removed later on
      return;

    int id = nuclib.NuclideIdsVsNames().at(mReactantNames[i]);
    chargeDiff -= mNsOfReactants[i] * nuclib.Zs()[id];
  }
  for (unsigned int i = 0; i < mProductNames.size(); ++i) {
    if (nuclib.NuclideIdsVsNames().count(mProductNames[i]) == 0)
      // ignore this reaction, it will be removed later on
      return;

    int id = nuclib.NuclideIdsVsNames().at(mProductNames[i]);
    chargeDiff += mNsOfProducts[i] * nuclib.Zs()[id];
  }

  if ((chargeDiff != 0) && !isWeak) {
    printf("Offending reaction (label %s): %s\n", mLabel.c_str(), ToString().c_str());
    throw std::invalid_argument("Attempted to create non-weak reaction with "
        "a charge difference in the nuclides");
  }

  if (chargeDiff == 0) {
    // there are no leptons
    return;
  }

  // always look at reaction going from Z to Z+1, if the charge difference is
  // negative, flip the lepton reactants and products and treat the reaction
  // as inverse
  bool inverse = isInverse;
  auto pLeptonReactants = &mReactantLeptons;
  auto pLeptonProducts = &mProductLeptons;

  if (chargeDiff < 0) {
    inverse = !isInverse;
    pLeptonReactants = &mProductLeptons;
    pLeptonProducts = &mReactantLeptons;
    chargeDiff = -chargeDiff;
  }

  for (int i = 0; i < chargeDiff; ++i) {
    if (inverse) {
      if (isDecay) {
        pLeptonReactants->push_back(Lepton::AntiE);
        pLeptonReactants->push_back(Lepton::NuE);
      } else {
        pLeptonReactants->push_back(Lepton::NuE);
        pLeptonProducts->push_back(Lepton::E);
      }
    } else {
      if (isDecay) {
        pLeptonProducts->push_back(Lepton::E);
        pLeptonProducts->push_back(Lepton::AntiNuE);
      } else {
        pLeptonReactants->push_back(Lepton::AntiE);
        pLeptonProducts->push_back(Lepton::AntiNuE);
      }
    }
  }

  mString = ToString();
}

Reaction::Reaction(const std::vector<std::string>& reactantNames,
    const std::vector<std::string>& productNames,
    const std::vector<int>& NsOfReactants,
    const std::vector<int>& NsOfProducts,
    const std::vector<Lepton>& reactantLeptons,
    const std::vector<Lepton>& productLeptons,
    const bool isWeak, const bool isInverse, const std::string& label) :
    mReactantNames(reactantNames),
    mProductNames(productNames),
    mNsOfReactants(NsOfReactants),
    mNsOfProducts(NsOfProducts),
    mReactantLeptons(reactantLeptons),
    mProductLeptons(productLeptons),
    mIsWeak(isWeak),
    mIsInverse(isInverse),
    mLabel(label) {
  Init();
}

void Reaction::AddUsedNuclideNames(
    std::unordered_set<std::string> * const pUsedNames) const {
  for (auto reactant : mReactantNames)
    pUsedNames->insert(reactant);
  for (auto product : mProductNames)
    pUsedNames->insert(product);
}

bool Reaction::ContainsOnlyKnownNuclides(const NuclideLibrary& nuclib) const {
  for (auto iso : mReactantNames) {
    auto it = std::find(nuclib.Names().begin(), nuclib.Names().end(), iso);
    if (it == nuclib.Names().end()) return false;
  }
  for (auto iso : mProductNames) {
    auto it = std::find(nuclib.Names().begin(), nuclib.Names().end(), iso);
    if (it == nuclib.Names().end()) return false;
  }
  return true;
}

Reaction Reaction::Inverse() const {
  return Reaction(mProductNames, mReactantNames, mNsOfProducts, mNsOfReactants,
      mProductLeptons, mReactantLeptons, mIsWeak, !mIsInverse, mLabel);
}

bool Reaction::Less(const Reaction& other) const {
  return mString < other.mString;
}

bool Reaction::IsEquivalent(const Reaction& other) const {
  if (mReactantNames != other.mReactantNames) return false;
  if (mProductNames != other.mProductNames) return false;
  if (mNsOfReactants != other.mNsOfReactants) return false;
  if (mNsOfProducts != other.mNsOfProducts) return false;
  if (mReactantLeptons != other.mReactantLeptons) return false;
  if (mProductLeptons != other.mProductLeptons) return false;
  if (mIsWeak != other.mIsWeak) return false;

  return true;
}

bool Reaction::IsDecay() const {
  if (!mIsWeak)
    return false;

  if (mIsInverse)
    return (mProductLeptons.size() == 0);
  else
    return (mReactantLeptons.size() == 0);
}

void Reaction::Init() {
  if (mReactantNames.size() <= 0)
    throw std::invalid_argument("Attempted to create reaction with no "
        "reactants");

  if (mProductNames.size() <= 0)
    throw std::invalid_argument("Attempted to create reaction with no "
        "products");

  if (mReactantNames.size() != mNsOfReactants.size())
    throw std::invalid_argument("Attempted to create reaction with different "
        "number of reactant names and numbers");

  if (mProductNames.size() != mNsOfProducts.size())
    throw std::invalid_argument("Attempted to create reaction with different "
        "number of product names and numbers");

  // sort names and numbers
  SortNamesAndNumbers(&mReactantNames, &mNsOfReactants);
  SortNamesAndNumbers(&mProductNames, &mNsOfProducts);

  mDensityExponent = -1;
  for (int n : mNsOfReactants)
    mDensityExponent += n;

  mString = ToString();
}

std::string Reaction::ToString() const {
  std::ostringstream str;
  for (unsigned int i = 0; i < mReactantNames.size(); ++i) {
    for (int j = 0; j < mNsOfReactants[i]; ++j) {
      str << mReactantNames[i];
      if ((i != mReactantNames.size() - 1) || (j < mNsOfReactants[i] - 1))
        str << " + ";
    }
  }

  for (unsigned int i = 0; i < mReactantLeptons.size(); ++i) {
    str << " + ";
    str << Leptons::ToString(mReactantLeptons[i]);
  }

  str << " -> ";

  for (unsigned int i = 0; i < mProductNames.size(); ++i) {
    for (int j = 0; j < mNsOfProducts[i]; ++j) {
      str << mProductNames[i];
      if ((i != mProductNames.size() - 1) || (j < mNsOfProducts[i] - 1))
        str << " + ";
    }
  }

  for (unsigned int i = 0; i < mProductLeptons.size(); ++i) {
    str << " + ";
    str << Leptons::ToString(mProductLeptons[i]);
  }

  return str.str();
}
