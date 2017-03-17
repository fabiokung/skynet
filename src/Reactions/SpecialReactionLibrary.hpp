/// \file SpecialReactionLibrary.hpp
/// \author lroberts
/// \since Mar 17, 2017
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_SPECIALREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_SPECIALREACTIONLIBRARY_HPP_

#include <functional> 
#include <vector>
#include <utility> 

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"
#include "EquationsOfState/EOS.hpp"

typedef std::function<double(double, ThermodynamicState)> RateFunc; 

class SpecialReactionLibrary: public ReactionLibraryBase {
public:
  SpecialReactionLibrary(
      const std::vector< std::pair<Reaction, RateFunc> >& reactions,
      const ReactionType reacType, 
      const std::string& description, 
      const std::string& source, 
      const NuclideLibrary& nucLib,
      const NetworkOptions& opts);

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new SpecialReactionLibrary(*this));
  }

  bool StateChanged(const ThermodynamicState /*thermoState*/,
      NetworkOutput * const /*pOutput*/) {
    // For now there is no state change, since the reactions
    // are assumed constant when they go off the edge of the
    // table.  At low temperature and density this is a good
    // assumption because capture rates go to zero and the
    // decay rates are the free space decay rates.  At high
    // density and temperature this is not so good, but if
    // we are out of the allowed upper limits of the tables
    // these reactions should definitely *not* be set to
    // zero.
    return false;
  }

  std::string Name() const {
    return "Flexible Reaction Rates Library";
  }

  const std::vector<double>& Rates() const {
    return mRates;
  }

  const std::vector<double>& InverseRates() const {
    return mInverseRates;
  }

  const std::vector<double>& HeatingRates() const {
      return mHeatingRates;
  }

  const std::vector<double>& HeatingInverseRates() const {
      return mHeatingInverseRates;
  }

protected:
  void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func);

  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap);

private:
  std::vector<double> mRates;
  std::vector<double> mInverseRates;
  std::vector<double> mHeatingRates;
  std::vector<double> mHeatingInverseRates;

  ReactionData<RateFunc> mRateFuncs;   // Reaction Q-values
};

#endif // SKYNET_REACTIONS_SPECIALREACTIONLIBRARY_HPP_
