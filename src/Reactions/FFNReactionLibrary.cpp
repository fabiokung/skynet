/// \file FFNReactionLibrary.cpp
/// \author lroberts
/// \since Jan 21, 2015
///
/// \brief
///
///

#include "Reactions/FFNReactionLibrary.hpp"

#include <cmath>

FFNReactionLibrary::FFNReactionLibrary(const std::string& ffnPath,
    const ReactionType type, const std::string& description,
    const NuclideLibrary& nuclib, const NetworkOptions& opts,
    const bool doScreening) :
    FFNReactionLibrary(FFN(ffnPath, nuclib), type, description, nuclib, opts,
        doScreening) {}

FFNReactionLibrary::FFNReactionLibrary(const FFN& FFNlib,
    const ReactionType type, const std::string& description,
    const NuclideLibrary& nucLib, const NetworkOptions& opts,
    const bool doScreening) :
    ReactionLibraryBase(type, description, FFNlib.getSource(),
        FFNlib.GetValidReactions(nucLib), nucLib, opts, doScreening) {
  // load the YeRho, T9 grid
  mT9Grid = FFNlib.T9Grid();
  mLogYeRhoGrid = FFNlib.LogYeRhoGrid();

  // load the tables point by point
  for (unsigned int i = 0; i < mT9Grid.size(); i++) {
    for (unsigned int j = 0; j < mLogYeRhoGrid.size(); j++) {
      mTable.push_back(ReactionData<double>(
          FFNlib.GetTableEntryAll(j, i, nucLib)));
    }
  }

  mRates = std::vector<double>(NumAllReactions());
}

void FFNReactionLibrary::DoLoopOverReactionData(
    const std::function<void(GeneralReactionData * const)>& func) {
  for (auto& tab : mTable)
    func(&tab);

  mRates = std::vector<double>(Reactions().ActiveData().size());
}

void FFNReactionLibrary::DoCalculateRates(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection) {
  if ((pZs != nullptr) && (pScreeningChemicalPotentialCorrection != nullptr))
    DoCalculateRatesScreening<true>(thermoState,
        partitionFunctionsWithoutSpinTerms, expArgumentCap, pZs,
        pScreeningChemicalPotentialCorrection);
  else
    DoCalculateRatesScreening<false>(thermoState,
        partitionFunctionsWithoutSpinTerms, expArgumentCap, nullptr, nullptr);
}

template <bool SCR>
void FFNReactionLibrary::DoCalculateRatesScreening(
    const ThermodynamicState thermoState,
    const std::vector<double>& /*partitionFunctionsWithoutSpinTerms*/,
    const double /*expArgumentCap*/, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection) {
  // find our place in the table and get the weights
  auto rhoParams = FindInterpolationParams(mLogYeRhoGrid,
      log10(thermoState.Ye() * thermoState.Rho()));
  auto T9Params = FindInterpolationParams(mT9Grid, thermoState.T9());

  // put the reactions in the rate array
  int idx11 = GetTableIdx(T9Params.idxLower, rhoParams.idxLower);
  int idx21 = GetTableIdx(T9Params.idxUpper, rhoParams.idxLower);
  int idx12 = GetTableIdx(T9Params.idxLower, rhoParams.idxUpper);
  int idx22 = GetTableIdx(T9Params.idxUpper, rhoParams.idxUpper);

  double w11 = T9Params.weightLower * rhoParams.weightLower;
  double w21 = T9Params.weightUpper * rhoParams.weightLower;
  double w12 = T9Params.weightLower * rhoParams.weightUpper;
  double w22 = T9Params.weightUpper * rhoParams.weightUpper;

  for (unsigned int i = 0; i < mRates.size(); i++) {
    mRates[i] = pow(10.0,
        mTable[idx11][i] * w11 + mTable[idx21][i] * w21
            + mTable[idx12][i] * w12 + mTable[idx22][i] * w22);

    if (SCR) {
      ReactionData_t data = Data()[i];

      ReactionIds8_t ids;
      ids.Uint64[0] = Ids0()[i].Uint64;
      ids.Uint64[1] = Ids1()[i].Uint64;

      ReactionNs8_t Ns;
      Ns.Uint32[0] = Ns0()[i].Uint32;
      Ns.Uint32[1] = Ns1()[i].Uint32;

      double sumMuReact = 0.0;
      int sumZReact = 0.0;

      for (int k = 0; k < data.NumReactants; ++k) {
        int n = Ns.Ns[k];
        int Z = (*pZs)[ids.Ids[k]];
        sumMuReact += (double)n * (*pScreeningChemicalPotentialCorrection)[Z];
        sumZReact += n * Z;
      }

      mRates[i] *= exp(
          sumMuReact - (*pScreeningChemicalPotentialCorrection)[sumZReact]);
    }
  }
}

FFNReactionLibrary::InterpolationParams
FFNReactionLibrary::FindInterpolationParams(const std::vector<double>& table,
    const double value) const {
  InterpolationParams params;

  if (value >= table[table.size() - 1]) {
    params.idxUpper = table.size() - 1;
    params.idxLower = params.idxUpper - 1;
    params.weightUpper = 1.0;
    params.weightLower = 0.0;
  } else if (value <= table[0]) {
    params.idxUpper = 1;
    params.idxLower = 0;
    params.weightUpper = 0.0;
    params.weightLower = 1.0;
  } else {
    for (unsigned int i = 1; i < table.size(); i++) {
      if (value <= table[i]) {
        params.idxUpper = i;
        break;
      }
    }

    params.idxLower = params.idxUpper - 1;
    params.weightUpper = (value - table[params.idxLower])
        / (table[params.idxUpper] - table[params.idxLower]);
    params.weightLower = 1.0 - params.weightUpper;
  }

  return params;
}

// explicit template instantiation
template void FFNReactionLibrary::DoCalculateRatesScreening<false>(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection);

template void FFNReactionLibrary::DoCalculateRatesScreening<true>(
    const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection);
