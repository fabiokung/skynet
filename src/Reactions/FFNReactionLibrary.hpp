/// \file FFNReactionLibrary.hpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_FFNREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_FFNREACTIONLIBRARY_HPP_

#include <fstream>
#include <vector>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/FFN.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

class FFNReactionLibrary: public ReactionLibraryBase {
public:
  void Dump(const std::string& fname) const {
    std::ofstream ofs(fname.c_str(),
        std::ostream::trunc | std::ostream::binary);
    boost::archive::binary_oarchive outArch(ofs);
    outArch << *this;
  }

  static FFNReactionLibrary ReadFromDisk(const std::string& fname,
      const NetworkOptions& opts) {
    FFNReactionLibrary lib(opts);
    std::ifstream ifs(fname.c_str(),
        std::ostream::in | std::ostream::binary);
    boost::archive::binary_iarchive inArch(ifs);
    inArch >> lib;
    return lib;
  }

  FFNReactionLibrary(const std::string& ffnPath, const ReactionType type,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, const bool doScreening);

  FFNReactionLibrary(const FFN& FFNLib, const ReactionType type,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, const bool doScreening);

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new FFNReactionLibrary(*this));
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
    return "FFN Type Rates Reaction Library";
  }

  const std::vector<double>& Rates() const {
    return mRates;
  }

  const std::vector<double>& InverseRates() const {
    return mInverseRates;
  }

protected:
  void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func);

  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

  template<bool SCR>
  void DoCalculateRatesScreening(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

private:
  FFNReactionLibrary(const NetworkOptions& opts):
    ReactionLibraryBase(opts) {}

  friend class boost::serialization::access;
  template <class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & boost::serialization::base_object<ReactionLibraryBase>(*this);
    ar & mRates;
    ar & mInverseRates;
    ar & mTable;
    ar & mT9Grid;
    ar & mLogYeRhoGrid;
  }

  struct InterpolationParams {
    unsigned int idxLower, idxUpper;
    double weightLower, weightUpper;
  };

  InterpolationParams FindInterpolationParams(const std::vector<double>& table,
      const double value) const;

  inline int GetTableIdx(const int iT9, const int iLogYeRho) const {
    return iT9 * mLogYeRhoGrid.size() + iLogYeRho;
  }

  std::size_t mNumDecayReactions;

  std::vector<double> mRates;
  std::vector<double> mInverseRates;
  std::vector<ReactionData<double>> mTable; // FFN Table with beta or capture
                                            // data depending on the reaction
  std::vector<double> mT9Grid;
  std::vector<double> mLogYeRhoGrid;
};

#endif // SKYNET_REACTIONS_FFNREACTIONLIBRARY_HPP_
