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
#include <stdexcept>
#include <string>
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

  // Relative accuracy the rate integrals over this distribution can be
  // required to reach. An analytic distribution supports whatever the
  // quadrature can deliver; a tabulated one is only defined to the accuracy of
  // its own interpolation, and demanding more makes the adaptive integrator
  // fail with a roundoff error instead of returning its (converged) result.
  virtual double IntegrationRelativeError() const {
    return 1.0e-12;
  }

  // When a rate integral fails to converge, the reaction library zeroes the
  // rate -- fine for an analytic distribution, whose integrals can genuinely
  // fail for extreme parameters. A tabulated distribution's integral fails only
  // when its grid is too coarse to resolve the spectrum, a setup error that
  // must not be silently zeroed, so it makes the failure fatal instead.
  virtual bool IntegrationFailureIsFatal() const {
    return false;
  }

  double LocalT9() const {
    return mLocalTMeV;
  }

  virtual void SetLocalT9(const double T9) {
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
