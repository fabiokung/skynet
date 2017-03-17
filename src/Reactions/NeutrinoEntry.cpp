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
#include <stdexcept>

#include "NuclearData/Nuclide.hpp"
#include "Utilities/H5Helper.hpp"
#include "Utilities/Constants.hpp"

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

NeutrinoEntry NeutrinoEntry::FromVacuumBetaDecay(const Reaction& reac, 
    double vacRate, const NuclideLibrary& nuclib, bool inverse) { 
  if (reac.ReactantNames().size() > 1 || reac.ProductNames().size() > 1
      || reac.NsOfReactants()[0] > 1 || reac.NsOfProducts()[0] > 1) { 
    // Shouldn't be using this kind of reaction here 
    throw std::domain_error(
      "Can't use vacuum beta decay neutrino rates with more than "
      "one product and reactant."); 
  } 
  const int idReactant = nuclib.NuclideIdsVsNames().at(reac.ReactantNames()[0]);
  const int idProduct = nuclib.NuclideIdsVsNames().at(reac.ProductNames()[0]);
  
  const double Q = nuclib.MassExcessesInMeV()[idReactant] 
      - nuclib.MassExcessesInMeV()[idProduct];
  const double m = Constants::ElectronMassInMeV;

  // This is the analytic expression for the beta decay phase space
  const double phaseSpace = (15 * pow(m, 4) * Q * log((sqrt(Q*Q - m*m) + Q)/m) 
    + sqrt(Q*Q - m*m) * ( 2 * pow(Q,4) - 9*m*m*Q*Q - 8*pow(m,4)))/60/pow(m,5);     
  
  const double rate_const = log(2.0) / (mK);  
  const double M = vacRate /(rate_const * phaseSpace);  
  
  bool isNue = true; 
  if (nuclib.Zs()[idReactant] < nuclib.Zs()[idProduct]) isNue = false;
  
  // Pass the Q value with the wrong sign for consistency with 
  // neutrino reaction library 
  if (inverse) { 
    return NeutrinoEntry({nuclib.As()[idProduct]}, {nuclib.Zs()[idProduct]}, 
        {nuclib.As()[idReactant]}, {nuclib.Zs()[idReactant]}, Q, M, 0.0, !isNue);
  } 
  return NeutrinoEntry({nuclib.As()[idReactant]}, {nuclib.Zs()[idReactant]}, 
      {nuclib.As()[idProduct]}, {nuclib.Zs()[idProduct]}, -Q, M, 0.0, isNue);
 
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
