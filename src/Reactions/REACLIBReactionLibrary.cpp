/// \file REACLIBReactionLibrary.cpp
/// \author jlippuner
/// \since Jul 24, 2014
///
/// \brief
///
///

#include "Reactions/REACLIBReactionLibrary.hpp"

#include <cmath>
#include <limits>

#include "Reactions/ReactionDataTypes.hpp"
#include "Utilities/Constants.hpp"

namespace { // unnamed so this function can only be used in this file

inline bool IsRateTemperatureDependent(
    const REACLIBCoefficients_t& REACLIBCoefficients) {
  for (unsigned int i = 1; i < REACLIBCoefficients.size(); ++i) {
    if (REACLIBCoefficients[i] != 0.0)
      return true;
  }

  return false;
}

inline bool IsRateConst(const REACLIBCoefficients_t& REACLIBCoefficients,
    const int densityExponent) {
  if (IsRateTemperatureDependent(REACLIBCoefficients))
    return false;

  if (densityExponent != 0)
    return false;

  return true;
}

std::vector<Reaction> SortAndCheckReactions(
    const std::vector<Reaction>& reactions,
    const std::vector<REACLIBCoefficients_t>& REACLIBCoefficients,
    const bool calculateStrongInverseRates, std::size_t * const pNumNonConst,
    std::array<ReactionData<double>, NumREACLIBCoefficients> * const
    pRateFittingCoefficients) {
  std::vector<Reaction> nonConstReactions;
  std::vector<std::array<double, NumREACLIBCoefficients>> nonConstCoeff;
  std::vector<Reaction> constReactions;
  std::vector<std::array<double, NumREACLIBCoefficients>> constCoeff;

  for (unsigned int i = 0; i < reactions.size(); ++i) {
    // if we calculate the strong inverse rates, we throw out the inverse
    // reactions
    if (calculateStrongInverseRates && reactions[i].IsInverse())
      continue;

    if (IsRateConst(REACLIBCoefficients[i], reactions[i].DensityExponent())) {
      constReactions.push_back(reactions[i]);
      constCoeff.push_back(REACLIBCoefficients[i]);
    } else {
      nonConstReactions.push_back(reactions[i]);
      nonConstCoeff.push_back(REACLIBCoefficients[i]);
    }
  }

  // TODO sort reactions according to number of reactants, products and ids of
  // reactants and products
  std::vector<Reaction> sortedReactions(nonConstReactions);
  sortedReactions.insert(sortedReactions.end(), constReactions.begin(),
      constReactions.end());

  if (pNumNonConst != nullptr)
    *pNumNonConst = nonConstReactions.size();

  if (pRateFittingCoefficients != nullptr) {
    for (unsigned int k = 0; k < pRateFittingCoefficients->size(); ++k) {
      std::vector<double> coeffs(sortedReactions.size());

      for (unsigned int i = 0; i < nonConstCoeff.size(); ++i)
        coeffs[i] = nonConstCoeff[i][k];

      for (unsigned int i = 0; i < constCoeff.size(); ++i)
        coeffs[i + nonConstReactions.size()] = constCoeff[i][k];

      (*pRateFittingCoefficients)[k] = coeffs;
    }
  }

  return sortedReactions;
}

} // namespace [unnamed]

REACLIBReactionLibrary::REACLIBReactionLibrary(const std::string& reaclibPath,
    const ReactionType type, const bool calclateInverses,
    const LeptonMode leptonMode, const std::string& description,
    const NuclideLibrary& nuclib, const NetworkOptions& opts,
    const bool doScreening, const bool freezeAtMinT9) :
    REACLIBReactionLibrary(REACLIB(reaclibPath), reaclibPath, type,
        calclateInverses, leptonMode, description, nuclib, opts, doScreening,
        freezeAtMinT9) {}

REACLIBReactionLibrary::REACLIBReactionLibrary(const REACLIB& reaclib,
    const std::string& source, const ReactionType type,
    const bool calclateInverses, const LeptonMode leptonMode,
    const std::string& description, const NuclideLibrary& nuclib,
    const NetworkOptions& opts, const bool doScreening,
    const bool freezeAtMinT9) :
    REACLIBReactionLibrary(type, calclateInverses, description,  source,
        reaclib.GetReactions(type, nuclib, leptonMode),
        reaclib.GetCoefficients(type, nuclib), nuclib, opts, doScreening,
        freezeAtMinT9) {}

REACLIBReactionLibrary::REACLIBReactionLibrary(const ReactionType type,
    const bool calclateInverses, const std::string& description,
    const std::string& source, const std::vector<Reaction>& reactions,
    const std::vector<REACLIBCoefficients_t>& REACLIBCoefficients,
    const NuclideLibrary& nuclib, const NetworkOptions& opts,
    const bool doScreening, const bool freezeAtMinT9) :
    ReactionLibraryBase(type, description, source,
        SortAndCheckReactions(reactions, REACLIBCoefficients,
            calclateInverses, nullptr, nullptr),
        nuclib, opts, doScreening),
    mCalcInverses(calclateInverses),
    mFreezeAtMinT9(freezeAtMinT9),
    mTempDepReactionsEnabled(true) {
  if ((type == ReactionType::Weak) && mCalcInverses)
    throw std::invalid_argument("Cannot calculate inverses for weak "
        "reactions.");

  // run this again to get the number of non-const reactions and the correct
  // ordering of the fitting coefficients (it doesn't work to pass in the
  // the pointers in the constructor)
  std::size_t numNonConst = 0;
  SortAndCheckReactions(reactions, REACLIBCoefficients, mCalcInverses,
      &numNonConst, &mRateFittingCoefficients);

  std::vector<bool> rateIsTemperatureDependent(NumAllReactions());
  mRates = std::vector<double>(NumAllReactions());
  std::vector<double> partfs(NumAllReactions());
  std::vector<bool> isInverseReaction(NumAllReactions());

  int numInverse = 0;
  if (mCalcInverses) {
    numInverse = NumAllReactions();

    // we don't treat constant rate differently if we need to calculate inverse
    // rates (there probably aren't any constant rates anyway)
    numNonConst = NumAllReactions();
  }

  mNumActiveNonConstReactions = numNonConst;

  mInverseRates = std::vector<double>(numInverse);
  std::vector<double> gammas(numInverse);
  std::vector<int> deltaNs(numInverse);
  std::vector<double> inverseRateSpinCorrections(numInverse);

  for (unsigned int i = 0; i < NumAllReactions(); ++i) {
    REACLIBCoefficients_t coeffs;
    for (unsigned int k = 0; k < NumREACLIBCoefficients; ++k)
      coeffs[k] = mRateFittingCoefficients[k].AllData()[i];

    rateIsTemperatureDependent[i] = IsRateTemperatureDependent(coeffs);

    if (mCalcInverses) {
      double gamma = 1.0;
      int deltaN = 0;
      double inverseRateSpinCorrection = 1.0;

      ReactionData_t data = Data()[i];

      ReactionNs8_t Ns;
      ReactionIds8_t ids;

      Ns.Uint32[0] = Ns0()[i].Uint32;
      Ns.Uint32[1] = Ns1()[i].Uint32;
      ids.Uint64[0] = Ids0()[i].Uint64;
      ids.Uint64[1] = Ids1()[i].Uint64;

      for (int k = 0; k < data.NumReactants + data.NumProducts; ++k) {
        int N = Ns.Ns[k];
        int id = ids.Ids[k];

        double mass = nuclib.MassesInMeV()[id];
        double spinTerm = nuclib.GroundStateSpinTerms()[id];
        double Nd = (double)N;

        if (k < data.NumReactants) {
          // this is a reactant
          gamma *= pow(mass, Nd);
          deltaN += N;
          inverseRateSpinCorrection *= pow(spinTerm, Nd);
        } else {
          // this is a product
          gamma /= pow(mass, Nd);
          deltaN -= N;
          inverseRateSpinCorrection /= pow(spinTerm, Nd);
        }
      }

      gammas[i] = pow(gamma, 1.5);
      deltaNs[i] = deltaN;
      inverseRateSpinCorrections[i] = inverseRateSpinCorrection;
    } else {
      isInverseReaction[i] = Reactions()[i].IsInverse();
    }
  }

  mRateIsTemperatureDependent = rateIsTemperatureDependent;
  mPartfs = partfs;
  mIsInverseReaction = isInverseReaction;
  mGammas = gammas;
  mDeltaNs = deltaNs;
  mInverseRateSpinCorrections = inverseRateSpinCorrections;

  CalculateConstantRates();
}

void REACLIBReactionLibrary::PrintAdditionalInfo(
    NetworkOutput * const pOutput) const {
  pOutput->Log("#   Inverse rates from detailed balance: %s\n",
      mCalcInverses ? "yes" : "no");
  std::string yes = "yes (at T9 = "
      + std::to_string(Options().MinT9ForStrongReactions) + ")";
  pOutput->Log("#   Freeze rates at min T9: %s\n",
      mFreezeAtMinT9 ? yes.c_str(): "no");
}

bool REACLIBReactionLibrary::StateChanged(const ThermodynamicState thermoState,
    NetworkOutput * const pOutput) {
  bool changed = false;

  if (mTempDepReactionsEnabled
      && (thermoState.T9() < Options().MinT9ForStrongReactions)) {
    mTempDepReactionsEnabled = false;
    if (pOutput != nullptr)
      pOutput->Log("# [%s \"%s\"]:\n# Turning OFF temperature dependent "
          "reactions (temperature below threshold)\n", Name().c_str(),
          Description().c_str());

    changed = true;
  } else if (!mTempDepReactionsEnabled
      && (thermoState.T9() >= Options().MinT9ForStrongReactions)) {
    mTempDepReactionsEnabled = true;
    if (pOutput != nullptr)
      pOutput->Log("# [%s \"%s\"]:\n# Turning ON temperature dependent "
          "reactions (temperature above threshold)\n", Name().c_str(),
          Description().c_str());

    changed = true;
  }

  std::vector<bool> activeFlags(NumAllReactions(), true);

  if (!mFreezeAtMinT9 && changed) {
    for (unsigned int i = 0; i < NumAllReactions(); ++i) {
      activeFlags[i] =
          (mTempDepReactionsEnabled || (!mRateIsTemperatureDependent[i]));
    }

    UpdateActive(activeFlags);
  }

  return changed;
}

void REACLIBReactionLibrary::DoLoopOverReactionData(
      const std::function<void(GeneralReactionData * const)>& func) {
  auto numActive = Reactions().ActiveData().size();
  mRates = std::vector<double>(numActive);

  for (unsigned int k = 0; k < NumREACLIBCoefficients; ++k)
    func(&mRateFittingCoefficients[k]);

  func(&mRateIsTemperatureDependent);
  func(&mPartfs);
  func(&mIsInverseReaction);

  if (mCalcInverses) {
    func(&mGammas);
    func(&mDeltaNs);
    func(&mInverseRateSpinCorrections);
    mInverseRates = std::vector<double>(numActive);
  }

  CalculateConstantRates();
}

void REACLIBReactionLibrary::DoCalculateRates(
    const ThermodynamicState inputThermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection) {
  auto thermoState = inputThermoState;
  if (mFreezeAtMinT9 && !mTempDepReactionsEnabled)
    thermoState.SetT9(Options().MinT9ForStrongReactions);
  if ((pZs == nullptr) || (pScreeningChemicalPotentialCorrection == nullptr))
    DoCalculateRatesScreening<false>(thermoState,
        partitionFunctionsWithoutSpinTerms, expArgumentCap, nullptr, nullptr);
  else
    DoCalculateRatesScreening<true>(thermoState,
        partitionFunctionsWithoutSpinTerms, expArgumentCap, pZs,
        pScreeningChemicalPotentialCorrection);
}

template <bool SCR>
void REACLIBReactionLibrary::DoCalculateRatesScreening(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection) {
  CalculateREACLIBRateExpArguments(thermoState.T9());

  if (mCalcInverses) {
    double T9inMeV = thermoState.T9() * Constants::BoltzmannConstantInMeVPerGK;
    double T9TimesInverseRateFactor =
        thermoState.T9() * Constants::InverseRateFactor;

    // we only need to calculate the non-constant rates, constant rates have
    // exactly one reactant and hence no screening correction
    for (unsigned int i = 0; i < mNumActiveNonConstReactions; ++i) {
      double forwardRateExpArgument = mRates[i];

      ReactionData_t data = Data()[i];
      int densityExponent = (int)data.DensityExponent;

      ReactionFactorials_t factorials = Factorials()[i];
      double factorialProduct =
          (double)factorials.NsFactorialProdcutOfReactants;

      ReactionIds8_t ids;
      ids.Uint64[0] = Ids0()[i].Uint64;
      ids.Uint64[1] = Ids1()[i].Uint64;

      ReactionNs8_t Ns;
      Ns.Uint32[0] = Ns0()[i].Uint32;
      Ns.Uint32[1] = Ns1()[i].Uint32;

      double sumMuReact = 0.0;
      int sumZReact = 0;
      double sumMuProd = 0.0;
      int sumZProd = 0;
      double inverseRatePartitionFunctionFactor = 1.0;

      for (int k = 0; k < data.NumReactants; ++k) {
        int id = ids.Ids[k];
        int n = Ns.Ns[k];
        inverseRatePartitionFunctionFactor *=
            pow(partitionFunctionsWithoutSpinTerms[id], (double)n);

        if (SCR) {
          int Z = (*pZs)[id];
          sumMuReact += (double)n * (*pScreeningChemicalPotentialCorrection)[Z];
          sumZReact += n * Z;
        }
      }

      for (int k = data.NumReactants; k < data.NumReactants + data.NumProducts;
          ++k) {
        int id = ids.Ids[k];
        int n = Ns.Ns[k];
        inverseRatePartitionFunctionFactor /=
            pow(partitionFunctionsWithoutSpinTerms[id], (double)n);

        if (SCR) {
          int Z = (*pZs)[id];
          sumMuProd += (double)n * (*pScreeningChemicalPotentialCorrection)[Z];
          sumZProd += n * Z;
        }
      }

      if (SCR) {
        forwardRateExpArgument +=
            sumMuReact - (*pScreeningChemicalPotentialCorrection)[sumZReact];
      }

      // if the reaction is purely strong, sumZProd == sumZReact, but there are
      // some multi-species weak reactions where this may not be the case
      double expDiff = -Qs()[i] / T9inMeV + sumMuProd - sumMuReact;
      if (SCR) {
        expDiff += -(*pScreeningChemicalPotentialCorrection)[sumZProd]
          + (*pScreeningChemicalPotentialCorrection)[sumZReact];
      }

      // cap the arguments of the forward and inverse rate so that they are
      // still consistent with each other
      if (forwardRateExpArgument > expArgumentCap)
        forwardRateExpArgument = expArgumentCap;

      double inverseRateExpArgument = forwardRateExpArgument + expDiff;

      if (inverseRateExpArgument > expArgumentCap) {
        inverseRateExpArgument = expArgumentCap;
        forwardRateExpArgument = inverseRateExpArgument - expDiff;
      }

      mRates[i] = exp(forwardRateExpArgument)
          * pow(thermoState.Rho(), (double)densityExponent) / factorialProduct;

      double inverseRate = exp(inverseRateExpArgument);
      inverseRate *= mGammas[i] * mInverseRateSpinCorrections[i];
      inverseRate *=
          pow(T9TimesInverseRateFactor, double(3 * mDeltaNs[i]) / 2.0);

      inverseRate *= pow(thermoState.Rho(),
          (double)(densityExponent - mDeltaNs[i]));
      inverseRate /= factorialProduct;

      inverseRate *= inverseRatePartitionFunctionFactor;

      mInverseRates[i] = inverseRate;
    }
  } else {
    // note that mNumActiveNonConstReactions <= NumActiveReactions()
    for (unsigned int i = 0; i < NumActiveReactions(); ++i) {
      double partitionFunction = 1.0;
      ReactionData_t data;
      ReactionIds8_t ids;
      ReactionNs8_t Ns;

      double sumMu = 0.0;
      int sumZ = 0.0;

      if ((SCR && (i < mNumActiveNonConstReactions)) || mIsInverseReaction[i]) {
        data = Data()[i];

        ids.Uint64[0] = Ids0()[i].Uint64;
        ids.Uint64[1] = Ids1()[i].Uint64;

        Ns.Uint32[0] = Ns0()[i].Uint32;
        Ns.Uint32[1] = Ns1()[i].Uint32;

        for (int k = 0; k < data.NumReactants; ++k) {
          int id = ids.Ids[k];
          int n = Ns.Ns[k];

          if (mIsInverseReaction[i]) {
            partitionFunction *=
                pow(partitionFunctionsWithoutSpinTerms[id], (double)n);
          }

          if (SCR) {
            int Z = (*pZs)[id];
            sumMu += (double)n * (*pScreeningChemicalPotentialCorrection)[Z];
            sumZ += n * Z;
          }
        }

        if (mIsInverseReaction[i]) {
          for (int k = data.NumReactants;
              k < data.NumReactants + data.NumProducts; ++k) {
            int id = ids.Ids[k];
            int n = Ns.Ns[k];
            partitionFunction /=
                pow(partitionFunctionsWithoutSpinTerms[id], (double)n);
          }

        }
      }

      if (i < mNumActiveNonConstReactions) {
        ReactionData_t data = Data()[i];
        double densityExponent = (double)data.DensityExponent;

        ReactionFactorials_t factorials = Factorials()[i];
        double factorialProduct =
            (double)factorials.NsFactorialProdcutOfReactants;

        double arg = mRates[i];

        if (SCR) {
          arg += sumMu - (*pScreeningChemicalPotentialCorrection)[sumZ];
        }

        if (arg > expArgumentCap) {
          mRates[i] = exp(expArgumentCap);
        } else {
          mRates[i] = exp(arg) * pow(thermoState.Rho(), densityExponent)
              / factorialProduct;
        }
      }

      if (mIsInverseReaction[i])
        // need to divide here because this is the inverse reaction, so
        // reactants and products are swapped
        mRates[i] /= partitionFunction;
    }
  }

  //DoCheckRates(mRates, false, warningThreshold);
  //DoCheckRates(mInverseRates, true, warningThreshold);
}

void REACLIBReactionLibrary::CalculateREACLIBRateExpArguments(const double T9) {
  double T9ToMinus1 = 1.0 / T9;
  double T9To1Third = pow(T9, 1.0 / 3.0);
  double T9ToMinus1Third = 1.0 / T9To1Third;
  double T9To5Thirds = T9 * T9To1Third * T9To1Third;
  double logT9 = log(T9);

  // we only need to calculate the non-constant rates
  for (unsigned int i = 0; i < mNumActiveNonConstReactions; ++i) {
    double arg = mRateFittingCoefficients[0][i];
    arg += mRateFittingCoefficients[1][i] * T9ToMinus1;
    arg += mRateFittingCoefficients[2][i] * T9ToMinus1Third;
    arg += mRateFittingCoefficients[3][i] * T9To1Third;
    arg += mRateFittingCoefficients[4][i] * T9;
    arg += mRateFittingCoefficients[5][i] * T9To5Thirds;
    arg += mRateFittingCoefficients[6][i] * logT9;

    mRates[i] = arg;
  }
}

void REACLIBReactionLibrary::CalculateConstantRates() {
  // first find number of non constant active reactions and make sure they are
  // all at the beginning
  if (!mCalcInverses) {
    std::size_t numNonConst = 0;
    bool foundConstReaction = false;

    for (unsigned int i = 0; i < NumActiveReactions(); ++i) {
      REACLIBCoefficients_t coeffs;
      for (unsigned int k = 0; k < NumREACLIBCoefficients; ++k)
        coeffs[k] = mRateFittingCoefficients[k][i];

      bool isConst = IsRateConst(coeffs, Reactions()[i].DensityExponent());

      if (isConst) {
        if (!foundConstReaction)
          foundConstReaction = true;
      } else {
        if (foundConstReaction)
          throw std::runtime_error("Found a non-const reaction after a const "
              "reaction: " + Reactions()[i].String());
        ++numNonConst;
      }
    }

    mNumActiveNonConstReactions = numNonConst;
  } else {
    mNumActiveNonConstReactions = NumActiveReactions();
  }

  for (unsigned int i = mNumActiveNonConstReactions; i < NumActiveReactions();
      ++i) {
    ReactionFactorials_t factorials = Factorials()[i];
    double factorialProduct =
        (double)factorials.NsFactorialProdcutOfReactants;
    mRates[i] = exp(mRateFittingCoefficients[0][i]) / factorialProduct;
  }
}

// explicit template instantiation
template void REACLIBReactionLibrary::DoCalculateRatesScreening<false>(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection);

template void REACLIBReactionLibrary::DoCalculateRatesScreening<true>(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection);
