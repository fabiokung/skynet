/// \file REACLIBEntry.hpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#ifndef SKYNET_NUCLEARDATA_REACLIBENTRY_HPP_
#define SKYNET_NUCLEARDATA_REACLIBENTRY_HPP_

#include <array>
#include <string>
#include <vector>

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name RateType in Python instead of the enum values being
// global constants

#ifdef SWIG
%rename(RateType) RateTypeStruct;
#endif // SWIG

struct RateTypeStruct {
  enum Value {
    NonResonant,
    Resonant,
    Weak,
    SpontaneousDecay
  };
};

typedef RateTypeStruct::Value RateType;

constexpr std::size_t NumREACLIBCoefficients = 7;
typedef std::array<double, NumREACLIBCoefficients> REACLIBCoefficients_t;

class REACLIBEntry {
public:
  REACLIBEntry(const std::vector<std::string>& reactants,
      const std::vector<std::string>& products, const std::string& label,
      const RateType type, const bool isInverseRate, const double q,
      const REACLIBCoefficients_t a);

  const REACLIBCoefficients_t& RateFittingCoefficients() const {
    return mRateFittingCoefficients;
  }

  bool IsInverseRate() const {
    return mIsInverseRate;
  }

  const std::vector<std::string>& Reactants() const {
    return mReactants;
  }

  const std::vector<std::string>& Products() const {
    return mProducts;
  }

  double Q() const {
    return mQ;
  }

  bool IsWeak() const {
    return (mType == RateType::Weak);
  }

  const std::string& Label() const {
    return mLabel;
  }

private:
  std::vector<std::string> mReactants;
  std::vector<std::string> mProducts;
  std::string mLabel;
  RateType mType;
  bool mIsInverseRate;
  double mQ;
  REACLIBCoefficients_t mRateFittingCoefficients;
};

#endif // SKYNET_NUCLEARDATA_REACLIBENTRY_HPP_
