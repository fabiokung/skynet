/// \file NeutrinoDistribution.hpp
/// \author lroberts
/// \since May 20, 2015
///
/// \brief
///
///

#ifndef SKYNET_EQUATIONSOFSTATE_NEUTRINODISTRIBUTION_HPP_
#define SKYNET_EQUATIONSOFSTATE_NEUTRINODISTRIBUTION_HPP_

#include <functional>
#include <memory>
#include <limits>
#include <valarray>
#include <vector>

#include "Utilities/Constants.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name NeutrinoSpecies in Python instead of the enum values
// being global constants

#ifdef SWIG
%rename(NeutrinoSpecies) NeutrinoSpeciesStruct;
#endif // SWIG

struct NeutrinoSpeciesStruct {
  enum Value {
    NuE,
    AntiNuE,
    NuX,
    NuMu,
    AntiNuMu,
    NuTau,
    AntiNuTau
  };
};

typedef NeutrinoSpeciesStruct::Value NeutrinoSpecies;

class NeutrinoDistribution {
public:
  NeutrinoDistribution(const std::valarray<NeutrinoSpecies>& species) :
      mSpecies(species),
      mLocalTMeV(0.0) {}

  virtual ~NeutrinoDistribution() {}

  static std::string NeutrinoSpeciesToStr(const NeutrinoSpecies spec) {
    switch (spec) {
    case NeutrinoSpecies::NuE:
      return "NuE";
    case NeutrinoSpecies::AntiNuE:
      return "AntiNuE";
    case NeutrinoSpecies::NuX:
      return "NuX";
    case NeutrinoSpecies::NuMu:
      return "NuMu";
    case NeutrinoSpecies::AntiNuMu:
      return "AntiNuMu";
    case NeutrinoSpecies::NuTau:
      return "NuTau";
    case NeutrinoSpecies::AntiNuTau:
      return "AntiNuTau";
    default:
      throw std::runtime_error("Unknown neutrino species");
    }
  }

  // returns a distribution function f(e), where e is the neutrino energy in MeV
  virtual std::function<double(double)> DistributionFunction(
      const NeutrinoSpecies species) const =0;

  double LocalT9() const {
    return mLocalTMeV;
  }

  void SetLocalT9(const double T9) {
    mLocalTMeV = T9 * Constants::BoltzmannConstantInMeVPerGK;
  }

  const std::valarray<NeutrinoSpecies>& Species() const {
    return mSpecies;
  }

protected:
  std::valarray<NeutrinoSpecies> mSpecies;
  double mLocalTMeV; // local fluid temperature in MeV
};


class DummyNeutrinoDistribution : public NeutrinoDistribution {
public:
  DummyNeutrinoDistribution() :
      NeutrinoDistribution({ }) {}

  std::function<double(double)> DistributionFunction(
      const NeutrinoSpecies /*species*/) const {
    return [] (const double /*enuInMeV*/) { return 0.0; };
  }
};


#endif // SKYNET_EQUATIONSOFSTATE_NEUTRINODISTRIBUTION_HPP_
