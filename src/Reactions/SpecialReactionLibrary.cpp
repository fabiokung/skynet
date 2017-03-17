/// \file SpecialReactionLibrary.cpp
/// \author lroberts
/// \since March 18, 2017
///
/// \brief
///
///

#include "Reactions/SpecialReactionLibrary.hpp"
#include "Utilities/Constants.hpp"
#include "Reactions/ReactionData.hpp"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>

SpecialReactionLibrary::SpecialReactionLibrary(
    const std::vector< std::pair<Reaction, RateFunc> >& reactions,
    const ReactionType reacType, 
    const std::string& description, 
    const std::string& source, 
    const NuclideLibrary& nucLib,
    const NetworkOptions& opts) :
    ReactionLibraryBase(reacType, description, source, 
        [&reactions]{
          std::vector<Reaction> reacs;
          for (auto& reac : reactions) reacs.push_back(reac.first);  
          return reacs;
        }(), nucLib, opts) {
  std::vector<RateFunc> rateFuncs; 
  
  for (auto& reac : reactions) rateFuncs.push_back(reac.second);  
  mRateFuncs = ReactionData<RateFunc>(rateFuncs); 

  mRates = std::vector<double>(NumAllReactions());
  mInverseRates = std::vector<double>();
  mHeatingRates = std::vector<double>(NumAllReactions());
  mHeatingInverseRates = std::vector<double>();
}

void SpecialReactionLibrary::DoLoopOverReactionData(
    const std::function<void(GeneralReactionData * const)>& func) {
  func(&mRateFuncs);
  mRates = std::vector<double>(Reactions().ActiveData().size());
  mInverseRates = std::vector<double>(Reactions().ActiveData().size());
  mHeatingRates = std::vector<double>(Reactions().ActiveData().size());
  mHeatingInverseRates = std::vector<double>(Reactions().ActiveData().size());
}

void SpecialReactionLibrary::DoCalculateRates(
    const ThermodynamicState thermoState,
    const std::vector<double>& /*partitionFunctionsWithoutSpinTerms*/,
    const double /*expArgumentCap*/) {
  
  double time = 0.0; 
  for (unsigned int i = 0; i < mRateFuncs.size(); ++i) { 
    mRates[i] = mRateFuncs[i](time, thermoState);
    //mInverseRates[i] = 0.0; 
    mHeatingRates[i] = 0.0; 
    //mHeatingInverseRates[i] = 0.0;  
  }
}

