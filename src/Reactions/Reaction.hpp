/// \file Reaction.hpp
/// \author jlippuner
/// \since Jun 30, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACTION_HPP_
#define SKYNET_REACTIONS_REACTION_HPP_

#include <unordered_set>
#include <valarray>
#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "NuclearData/NuclideLibrary.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name Lepton in Python instead of the enum values
// being global constants

#ifdef SWIG
%rename(Lepton) LeptonStruct;
#endif // SWIG

struct LeptonStruct {
  enum Value {
    E,
    AntiE,
    Mu,
    AntiMu,
    Tau,
    AntiTau,
    NuE,
    AntiNuE,
    NuX,
    NuMu,
    AntiNuMu,
    NuTau,
    AntiNuTau
  };
};

typedef LeptonStruct::Value Lepton;

struct Leptons {
  static std::string ToString(const Lepton lepton);

  // (e, Le, Lmu, Ltau) where e is electric charge, Le is electronic lepton
  // number, Lmu is muonic lepton number, and Ltau is tauonic lepton number
  static std::valarray<int> Charge(const Lepton lepton);
};

class Reaction {
public:
  Reaction() {}

  Reaction(const std::vector<std::string>& reactantNames,
      const std::vector<std::string>& productNames,
      const std::vector<int>& NsOfReactants,
      const std::vector<int>& NsOfProducts, const bool isWeak,
      const bool isInverse, const bool isDecay, const std::string& label,
      const NuclideLibrary& nuclib);

  Reaction(const std::vector<std::string>& reactantNames,
      const std::vector<std::string>& productNames,
      const std::vector<int>& NsOfReactants,
      const std::vector<int>& NsOfProducts,
      const std::vector<Lepton>& reactantLeptons,
      const std::vector<Lepton>& productLeptons,
      const bool isWeak, const bool isInverse, const std::string& label);

  void AddUsedNuclideNames(
      std::unordered_set<std::string> * const pUsedNames) const;

  bool ContainsOnlyKnownNuclides(const NuclideLibrary& nuclib) const;

  Reaction Inverse() const;

  bool Less(const Reaction& other) const;

  bool IsEquivalent(const Reaction& other) const;

  const std::vector<std::string>& ReactantNames() const {
    return mReactantNames;
  }

  const std::vector<std::string>& ProductNames() const {
    return mProductNames;
  }

  const std::vector<int>& NsOfReactants() const {
    return mNsOfReactants;
  }

  const std::vector<int>& NsOfProducts() const {
    return mNsOfProducts;
  }

  const std::vector<Lepton> ReactantLeptons() const {
    return mReactantLeptons;
  }

  const std::vector<Lepton> ProductLeptons() const {
    return mProductLeptons;
  }

  int DensityExponent() const {
    return mDensityExponent;
  }

  bool IsWeak() const {
    return mIsWeak;
  }

  bool IsInverse() const {
    return mIsInverse;
  }

  bool IsDecay() const;

  const std::string& Label() const {
    return mLabel;
  }

  const std::string& String() const {
    return mString;
  }

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & mReactantNames;
    ar & mProductNames;
    ar & mNsOfReactants;
    ar & mNsOfProducts;
    ar & mReactantLeptons;
    ar & mProductLeptons;
    ar & mDensityExponent;
    ar & mIsWeak;
    ar & mIsInverse;
    ar & mLabel;
    ar & mString;
  }

  void Init();

  std::string ToString() const;

  // list of nuclide indices of the reactants and products
  std::vector<std::string> mReactantNames;
  std::vector<std::string> mProductNames;

  // the N_i^k integers of the reactants and products
  std::vector<int> mNsOfReactants;
  std::vector<int> mNsOfProducts;

  // list of leptons (possibly repeated if more than 1) participating in the
  // reaction
  std::vector<Lepton> mReactantLeptons;
  std::vector<Lepton> mProductLeptons;

  int mDensityExponent;

  bool mIsWeak;
  bool mIsInverse;
  std::string mLabel;

  std::string mString;
};

#ifndef SWIG
namespace std {
template<>
struct less<Reaction> {
  bool operator()(const Reaction& x, const Reaction& y) const {
    return x.Less(y);
  }
  typedef Reaction first_argument_type;
  typedef Reaction second_argument_type;
  typedef bool result_type;
};
} // namespace std
#endif // not SWIG

#endif // SKYNET_REACTIONS_REACTION_HPP_
