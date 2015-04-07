/// \file ReactionDataTypes.hpp
/// \author jlippuner
/// \since Jul 18, 2014
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACTIONDATATYPES_HPP_
#define SKYNET_REACTIONS_REACTIONDATATYPES_HPP_

#include <cstdint>

#include <boost/serialization/access.hpp>

#define CHECK_SIZE(type, size)                                                  \
static_assert(sizeof(type) == size, #type " has wrong size.");

struct ReactionData_t {
  // number of reactants
  uint8_t NumReactants;

  // number of products
  uint8_t NumProducts;

  // density exponent (sum_i N_i) - 1 where i ranges over reactants
  uint8_t DensityExponent;

  // (sum_i N_i) - (sum_j N_j), where i and j range over reactants and products,
  // respectively
  int8_t DeltaN;

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & NumReactants;
    ar & NumProducts;
    ar & DensityExponent;
    ar & DeltaN;
  }
};
CHECK_SIZE(ReactionData_t, 4);

union ReactionNs4_t {
  // This struct stores the numbers of the reactants and products that are
  // involved in a reaction for 4 different reactants or products.
  uint8_t Ns[4];
  uint32_t Uint32;

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & Uint32;
  }
};
CHECK_SIZE(ReactionNs4_t, 4);
BOOST_CLASS_IMPLEMENTATION(ReactionNs4_t, boost::serialization::object_serializable)

union ReactionNs8_t {
  uint8_t Ns[8];
  uint32_t Uint32[2];

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & Uint32;
  }
};
CHECK_SIZE(ReactionNs8_t, 8);
BOOST_CLASS_IMPLEMENTATION(ReactionNs8_t, boost::serialization::object_serializable)

//BOOST_CLASS_IMPLEMENTATION(uint64_t, boost::serialization::object_serializable)
union ReactionIds4_t {
  // This struct stores the ids of the reactants and products that are involved
  // in a reaction for 4 different reactants or products.
  uint16_t Ids[4];
  uint64_t Uint64;

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & Uint64;
  }
};
CHECK_SIZE(ReactionIds4_t, 8);
BOOST_CLASS_IMPLEMENTATION(ReactionIds4_t, boost::serialization::object_serializable)


union ReactionIds8_t {
  uint16_t Ids[8];
  uint64_t Uint64[2];

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & Uint64;
  }
};
CHECK_SIZE(ReactionIds8_t, 16);
BOOST_CLASS_IMPLEMENTATION(ReactionIds8_t, boost::serialization::object_serializable)

struct ReactionFactorials_t {
  // product of factorials of Ns of reactants
  uint16_t NsFactorialProdcutOfReactants;

  // product of factorials of Ns of products
  uint16_t NsFactorialProdcutOfProducts;

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & NsFactorialProdcutOfReactants;
    ar & NsFactorialProdcutOfProducts;
  }
};
CHECK_SIZE(ReactionFactorials_t, 4);

#endif // SKYNET_REACTIONS_REACTIONDATATYPES_HPP_
