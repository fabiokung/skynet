/// \file FFNReactionLibrary.hpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_NEUTRINOREACTIONLIBRARY_HPP_
#define SKYNET_REACTIONS_NEUTRINOREACTIONLIBRARY_HPP_

#include <vector>

#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Neutrino.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionLibraryBase.hpp"

// we use the workaround http://stackoverflow.com/a/13406244 to make the enum
// show up under the name NeutrinoCorrectionMode in Python instead of
// the enum values being global constants

#ifdef SWIG
%rename(NeutrinoCorrectionMode) NeutrinoCorrectionModeStruct;
#endif // SWIG

struct NeutrinoCorrectionModeStruct {
  enum Value {
    None,
    WeakMagnetism,
    WeakMagnetismExact
  };
};

typedef NeutrinoCorrectionModeStruct::Value NeutrinoCorrectionMode;

class NeutrinoReactionLibrary: public ReactionLibraryBase {
public:
  NeutrinoReactionLibrary(const Neutrino neutrinoLib,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, bool onlyNuCap = false,
      bool nuHeating = false, bool includeBeta = false,
      NeutrinoCorrectionMode correctionMode = NeutrinoCorrectionMode::None);

  NeutrinoReactionLibrary(const std::string& nuFile,
      const std::string& description, const NuclideLibrary& nuclib,
      const NetworkOptions& opts, bool onlyNuCap = false,
      bool nuHeating = false, bool includeBeta = false,
      NeutrinoCorrectionMode correctionMode = NeutrinoCorrectionMode::None) :
      NeutrinoReactionLibrary(Neutrino(nuFile, nuclib), description, nuclib,
        opts, onlyNuCap, nuHeating, includeBeta, correctionMode) {}

  std::unique_ptr<ReactionLibraryBase> MakeUniquePtr() const {
    return std::unique_ptr<ReactionLibraryBase>(
        new NeutrinoReactionLibrary(*this));
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
    return "Neutrino Reaction Rates Library";
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

  virtual void PrintAdditionalInfo(NetworkOutput * const pOutput) const;

protected:
  void DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func);

  void DoCalculateRates(const ThermodynamicState thermoState,
      const std::vector<double>& partitionFunctionsWithoutSpinTerms,
      const double expArgumentCap, const std::vector<int> * const pZs,
      const std::vector<double> * const pScreeningChemicalPotentialCorrection);

private:
  std::vector<double> mRates;
  std::vector<double> mInverseRates;
  std::vector<double> mHeatingRates;
  std::vector<double> mHeatingInverseRates;

  // Q values are defined in the direction of neutrino capture
  // so that E_nu = E_e - Q in both directions
  // e + A -> B + nu 
  // E_e + M_A = M_B + E_nu -> Q = E_e - E_nu = M_B - M_A 
  // which is opposite of the normal way of defining the Q-value 
  // considering that we are treating e + A -> B + nu as the forward reaction
  ReactionData<double> mQ;   // Reaction Q-values
  ReactionData<double> mMatrixElement; // Reaction matrix elements
  ReactionData<double> mWm;  // "Weak Magnetism" correction energy scale
  ReactionData<bool> mNueCap; // Is this an electron neutrino capture reaction
  const double mK = 6144.0; // Beta decay timescale constant (see Arcones 2010)
  bool mOnlyNuCap;
  bool mNuHeating;
  bool mIncludeBeta;
  NeutrinoCorrectionMode mCorrectionMode;
};

#endif // SKYNET_REACTIONS_NEUTRINOREACTIONLIBRARY_HPP_
