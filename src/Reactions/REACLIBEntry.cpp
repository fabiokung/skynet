/// \file Entry.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include "Reactions/REACLIBEntry.hpp"

#include <stdexcept>

REACLIBEntry::REACLIBEntry(const std::vector<std::string>& reactants,
    const std::vector<std::string>& products, const std::string& label,
    const RateType type, const bool isInverseRate, const double q,
    const REACLIBCoefficients_t a) :
    mReactants(reactants),
    mProducts(products),
    mLabel(label),
    mType(type),
    mIsInverseRate(isInverseRate),
    mQ(q),
    mRateFittingCoefficients(a) { }
