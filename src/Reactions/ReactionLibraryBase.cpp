/// \file ReactionLibraryBase.cpp
/// \author jlippuner
/// \since Jul 18, 2014
///
/// \brief
///
///

#include "Reactions/ReactionLibraryBase.hpp"

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Utilities/File.hpp"

namespace { // unnamed so these functions can only be used in this file

int Factorial(const int n) {
  assert(n >= 0);

  int factorial = 1;
  for (int i = 2; i <= n; ++i)
    factorial *= i;

  return factorial;
}

} // namespace [unnamed]

ReactionLibraryBase::ReactionLibraryBase(const ReactionType type,
    const std::string& description, const std::string& source,
    const std::vector<Reaction>& reactions, const NuclideLibrary& nuclib,
    const NetworkOptions& opts, const bool doScreening) :
    mType(type),
    mDescription(description),
    mSource(File::MakeAbsolutePath(source)),
    mStateHasBeenSet(false),
    mDoScreening(doScreening),
    mReactions(reactions),
    mDoScreenOutput(false),
    mOpts(opts),
    mEmptyDouble(std::vector<double>(0)) {
  std::vector<ReactionData_t> data(reactions.size());
  std::vector<ReactionNs4_t> ns0(reactions.size());
  std::vector<ReactionNs4_t> ns1(reactions.size());
  std::vector<ReactionIds4_t> ids0(reactions.size());
  std::vector<ReactionIds4_t> ids1(reactions.size());
  std::vector<ReactionFactorials_t> factorials(reactions.size());
  std::vector<double> qs(reactions.size());
  std::vector<std::string> labels(reactions.size());

  for (unsigned int i = 0; i < reactions.size(); ++i) {
    const Reaction& reaction = reactions[i];

    switch (mType) {
    case ReactionType::Weak:
      if (!reaction.IsWeak())
        throw std::invalid_argument("Passed non-weak reaction to a weak "
            "reaction library.");
      break;
    case ReactionType::Strong:
      if (reaction.IsWeak())
        throw std::invalid_argument("Passed weak reaction to a strong "
            "reaction library.");
      break;
    default:
      throw std::runtime_error("Unknown reaction type");
    }

    // add ReactionData_t and ReactionFactorials_t
    ReactionData_t reactionData;
    ReactionFactorials_t reactionFactorials;

    reactionData.NumReactants = reaction.ReactantNames().size();
    reactionData.NumProducts = reaction.ProductNames().size();

    if ((reactionData.NumReactants + reactionData.NumProducts) > 8)
      throw std::invalid_argument("Tried to add reaction '"
          + reaction.String() + "' with "
          + std::to_string(reactionData.NumReactants + reactionData.NumProducts)
          + " reactants and products to ReactionLibrary");

    reactionData.DensityExponent = 0;

    reactionFactorials.NsFactorialProdcutOfReactants = 1;
    reactionFactorials.NsFactorialProdcutOfProducts = 1;

    reactionData.DeltaN = 0;
    for (unsigned int k = 0; k < reaction.NsOfReactants().size(); ++k) {
      int N = reaction.NsOfReactants()[k];
      if (N < 1)
        throw std::invalid_argument("Reaction '" + reaction.String()
            + "' has a negative number of reactants");

      reactionData.DensityExponent += N;
      reactionData.DeltaN += N;
      reactionFactorials.NsFactorialProdcutOfReactants *= Factorial(N);
    }
    for (unsigned int k = 0; k < reaction.NsOfProducts().size(); ++k) {
      int N = reaction.NsOfProducts()[k];
      if (N < 1)
        throw std::invalid_argument("Reaction '" + reaction.String()
            + "' has a negative number of products");

      reactionData.DeltaN -= N;
      reactionFactorials.NsFactorialProdcutOfProducts *= Factorial(N);
    }

    reactionData.DensityExponent -= 1;

    data[i] = reactionData;
    factorials[i] = reactionFactorials;

    // add ReactionNs4_t and ReactionIds4_t
    ReactionNs8_t Ns;
    ReactionIds8_t ids;

    double Q = 0.0;
    int deltaA = 0;
    std::valarray<int> deltaCharge{0, 0, 0, 0};

    for (int k = 0; k < reactionData.NumReactants; ++k) {
      int N = reaction.NsOfReactants()[k];
      int id = nuclib.NuclideIdsVsNames().at(reaction.ReactantNames()[k]);
      double mass = nuclib.MassesInMeV()[id];
      int A = nuclib.As()[id];
      int Z = nuclib.Zs()[id];

      Ns.Ns[k] = N;
      ids.Ids[k] = id;
      Q += (double)N * mass;
      deltaA -= N * A;
      deltaCharge[0] += N * Z;
    }

    for (int k = 0; k < reactionData.NumProducts; ++k) {
      int N = reaction.NsOfProducts()[k];
      int id = nuclib.NuclideIdsVsNames().at(reaction.ProductNames()[k]);
      double mass = nuclib.MassesInMeV()[id];
      int A = nuclib.As()[id];
      int Z = nuclib.Zs()[id];

      Ns.Ns[reactionData.NumReactants + k] = N;
      ids.Ids[reactionData.NumReactants + k] = id;
      Q -= (double)N * mass;
      deltaA += N * A;
      deltaCharge[0] -= N * Z;
    }

    for (unsigned int k = 0; k < reaction.ReactantLeptons().size(); ++k) {
      deltaCharge += Leptons::Charge(reaction.ReactantLeptons()[k]);
    }

    for (unsigned int k = 0; k < reaction.ProductLeptons().size(); ++k) {
      deltaCharge -= Leptons::Charge(reaction.ProductLeptons()[k]);
    }

    if (deltaA != 0)
      throw std::invalid_argument("Reaction " + reaction.String() +
          " does not conserve mass");

    for (unsigned int k = 0; k < deltaCharge.size(); ++k) {
      if (deltaCharge[k] != 0) {
        throw std::invalid_argument("Reaction " + reaction.String() +
            " does not conserve charge or lepton number");
      }
    }

    ns0[i].Uint32 = Ns.Uint32[0];
    ns1[i].Uint32 = Ns.Uint32[1];
    ids0[i].Uint64 = ids.Uint64[0];
    ids1[i].Uint64 = ids.Uint64[1];

    qs[i] = Q;

    labels[i] = reaction.Label();
  }

  // make reaction data
  mData = data;
  mActiveFlags = std::vector<bool>(reactions.size(), true),
  mNs0 = ns0;
  mNs1 = ns1;
  mIds0 = ids0;
  mIds1 = ids1;
  mFactorials = factorials;
  mQs = qs;
  mLabels = labels;
}

void ReactionLibraryBase::CalculateRates(const ThermodynamicState thermoState,
    const std::vector<double>& partitionFunctionsWithoutSpinTerms,
    const double expArgumentCap, const std::vector<int> * const pZs,
    const std::vector<double> * const pScreeningChemicalPotentialCorrection) {
  if ((pZs == nullptr) != (pScreeningChemicalPotentialCorrection == nullptr))
    throw std::runtime_error("Zs and ScreeningChemicalPotentialCorrection must "
        "either both be given or both be nullptr.");

  if (!mStateHasBeenSet)
    SetState(thermoState);

  if (mDoScreening)
    DoCalculateRates(thermoState, partitionFunctionsWithoutSpinTerms,
        expArgumentCap, pZs, pScreeningChemicalPotentialCorrection);
  else
    DoCalculateRates(thermoState, partitionFunctionsWithoutSpinTerms,
        expArgumentCap, nullptr, nullptr);
}

void ReactionLibraryBase::RemoveReactions(
    const std::vector<Reaction>& reactionsToBeRemoved) {
  std::vector<bool> removeFlags(NumAllReactions(), false);

  // set vector of reactions to remove
  for (unsigned int i = 0; i < NumAllReactions(); i++) {
    removeFlags[i] = std::any_of(reactionsToBeRemoved.begin(),
        reactionsToBeRemoved.end(), [&] (const Reaction reac) {
          return reac.IsEquivalent(mReactions[i]);
        });
    // if (!remove[i]) std::cout << mReactions[i].GetReactionName() << "\n";
  }
  //std::cout << "Removing " << std::accumulate(remove.begin(),remove.end(),0,[](int a,bool b){if (!b) a++; return a;}) << " reactions. \n";

  // remove data associated with these reactions
  LoopOverReactionDataInBase([&] (GeneralReactionData * const pData) {
    pData->RemoveData(removeFlags);
  });

  DoLoopOverReactionData([&] (GeneralReactionData * const pData) {
    pData->RemoveData(removeFlags);
  });
}

void ReactionLibraryBase::RemoveReactions(
    const ReactionData<Reaction>& reactionsToBeRemoved) {
  RemoveReactions(reactionsToBeRemoved.AllData());
}

void ReactionLibraryBase::SetState(const ThermodynamicState thermoState) {
  mStateHasBeenSet = true;
  StateChanged(thermoState, nullptr);
}

void ReactionLibraryBase::ReIndexNuclides(const NuclideLibrary& nuclib) {
  std::vector<ReactionIds4_t> ids0(NumAllReactions());
  std::vector<ReactionIds4_t> ids1(NumAllReactions());

  for (unsigned int i = 0; i < NumAllReactions(); ++i) {
    ReactionIds8_t ids;
    auto data = mData[i];

    for (int k = 0; k < data.NumReactants; ++k)
      ids.Ids[k] =
          nuclib.NuclideIdsVsNames().at(mReactions[i].ReactantNames()[k]);

    for (int k = 0; k < data.NumProducts; ++k)
      ids.Ids[data.NumReactants + k] =
          nuclib.NuclideIdsVsNames().at(mReactions[i].ProductNames()[k]);

    ids0[i].Uint64 = ids.Uint64[0];
    ids1[i].Uint64 = ids.Uint64[1];
  }

  mIds0 = ids0;
  mIds1 = ids1;

  mIds0.SetActive(mActiveFlags.AllData());
  mIds1.SetActive(mActiveFlags.AllData());
}

void ReactionLibraryBase::AddYdotContributions(const std::vector<double>& Y,
    std::vector<double> * const pYdot) const {
  DoAddYdotContributions(Y, pYdot, Rates(), InverseRates());
}

double ReactionLibraryBase::GetEDotExternal(const std::vector<double>& Y) const {
  return DoGetEDotExternal(Y, HeatingRates(), HeatingInverseRates());
}

void ReactionLibraryBase::AddIndicesOfNonZeroJacobianEntries(
    std::set<std::array<int, 2>> * const indices) const {
  ReactionLibraryBase::DoAddIndicesOfNonZeroJacobianEntries(indices,
      (InverseRates().size() == Rates().size()));
}

void ReactionLibraryBase::AddJacobianContributions(
    const std::vector<double>& Y, Jacobian * const pJ) const {
  ReactionLibraryBase::DoAddJacobianContributions(Y, pJ, Rates(),
      InverseRates());
}

void ReactionLibraryBase::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# %s\n", Name().c_str());

  std::string typePrefix = "#   Reaction type: ";
  switch (mType) {
  case ReactionType::Weak:
    pOutput->Log("%sWeak\n", typePrefix.c_str());
    break;
  case ReactionType::Strong:
    pOutput->Log("%sStrong\n", typePrefix.c_str());
    break;
  default:
    throw std::runtime_error("Unknown reaction type");
  }

  pOutput->Log("#   Description: %s\n", Description().c_str());
  pOutput->Log("#   Source: %s\n", Source().c_str());
  pOutput->Log("#   Screening: %s\n", mDoScreening ? "yes" : "no");
  pOutput->Log("#   Forward reactions: %lu\n", Rates().size());
  pOutput->Log("#   Inverse reactions: %lu\n", InverseRates().size());
  PrintAdditionalInfo(pOutput);
  pOutput->Log("#\n");
}

void ReactionLibraryBase::PrintRates(NetworkOutput * const pOutput) const {
  PrintInfo(pOutput);

  pOutput->Log("#   Forward rates (may contain inverses if they are not "
      "specifically\n"
      "#                  calculated as inverse rates):\n");
  for (unsigned int i = 0; i < NumActiveReactions(); ++i) {
    std::string desc = mReactions[i].String();
    pOutput->Log("#   %25s: % 25.15E\n", desc.c_str(), Rates()[i]);
  }

  if (InverseRates().size() > 0) {
    pOutput->Log("#   Inverse rates:\n");
    for (unsigned int i = 0; i < NumActiveReactions(); ++i) {
      std::string desc = mReactions[i].Inverse().String();
      pOutput->Log("#   %25s: % 25.15E\n", desc.c_str(), InverseRates()[i]);
    }
  }
}

int ReactionLibraryBase::MaxScreeningZ(const std::vector<int>& Zs) {
  if (!mDoScreening)
    return 0;

  int maxZ = 0;

  for (unsigned int k = 0; k < NumActiveReactions(); ++k) {
    ReactionData_t data = mData[k];
    ReactionIds8_t ids;
    ids.Uint64[0] = mIds0[k].Uint64;
    ids.Uint64[1] = mIds1[k].Uint64;
    ReactionNs8_t Ns;
    Ns.Uint32[0] = mNs0[k].Uint32;
    Ns.Uint32[1] = mNs1[k].Uint32;

    int sumZ = 0.0;
    for (int m = 0; m < data.NumReactants; ++m) {
      int Z = Zs[ids.Ids[m]];
      maxZ = std::max(maxZ, Z);
      sumZ += Ns.Ns[m] * Z;
    }
    maxZ = std::max(maxZ, sumZ);

    sumZ = 0;
    for (int m = data.NumReactants; m < data.NumReactants + data.NumProducts;
        ++m) {
      int Z = Zs[ids.Ids[m]];
      maxZ = std::max(maxZ, Z);
      sumZ += Ns.Ns[m] * Z;
    }
    maxZ = std::max(maxZ, sumZ);
  }

  return maxZ;
}

void ReactionLibraryBase::UpdateActive(const std::vector<bool>& activeFlags) {
  if (activeFlags.size() != NumAllReactions())
    throw std::invalid_argument("Called UpdateActive with wrong number of "
        "flags");

  mActiveFlags = activeFlags;

  LoopOverReactionDataInBase([&] (GeneralReactionData * const pData) {
    pData->SetActive(mActiveFlags.AllData());
  });

  DoLoopOverReactionData([&] (GeneralReactionData * const pData) {
    pData->SetActive(mActiveFlags.AllData());
  });
}

void ReactionLibraryBase::LoopOverReactionDataInBase(
    const std::function<void(GeneralReactionData * const)>& func) {
  func(&mReactions);
  func(&mActiveFlags);
  func(&mData);
  func(&mNs0);
  func(&mNs1);
  func(&mIds0);
  func(&mIds1);
  func(&mFactorials);
  func(&mQs);
  func(&mLabels);
}

double ReactionLibraryBase::DoGetEDotExternal(const std::vector<double>& Y,
    const std::vector<double>& heatingRates,
    const std::vector<double>& heatingInverseRates) const {
  if (heatingRates.size() + heatingInverseRates.size() == 0) return 0.0;

  double edot = 0.0;
  bool alsoDoInverse = heatingInverseRates.size() > 0;
  for (unsigned int k = 0; k < NumActiveReactions(); ++k) {
    ReactionData_t data = mData[k];
    ReactionIds8_t ids;
    ids.Uint64[0] = mIds0[k].Uint64;
    ids.Uint64[1] = mIds1[k].Uint64;
    ReactionNs8_t Ns;
    Ns.Uint32[0] = mNs0[k].Uint32;
    Ns.Uint32[1] = mNs1[k].Uint32;

    double factorForward = heatingRates[k];
    double factorInverse = 0.0;

    for (int m = 0; m < data.NumReactants; ++m) {
      factorForward *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
    }

    if (alsoDoInverse) {
      factorInverse = heatingInverseRates[k];
      for (int m = data.NumReactants; m < data.NumReactants + data.NumProducts;
          ++m) {
        factorInverse *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
      }
    }

    edot += factorForward + factorInverse;
  }

  return edot;

}

void ReactionLibraryBase::DoAddYdotContributions(const std::vector<double>& Y,
    std::vector<double> * const pYdot, const std::vector<double>& rates,
    const std::vector<double>& inverseRates) const {
  bool alsoDoInverse = inverseRates.size() > 0;

  for (unsigned int k = 0; k < NumActiveReactions(); ++k) {
    ReactionData_t data = mData[k];
    ReactionIds8_t ids;
    ids.Uint64[0] = mIds0[k].Uint64;
    ids.Uint64[1] = mIds1[k].Uint64;
    ReactionNs8_t Ns;
    Ns.Uint32[0] = mNs0[k].Uint32;
    Ns.Uint32[1] = mNs1[k].Uint32;

    // these factors include the rate and the product of Y_m^N_m, where m ranges
    // over the reactants (or the product for the inverse)
    double factorForward = rates[k];
    double factorInverse = 0.0;
    if (alsoDoInverse)
      factorInverse = inverseRates[k];

    // m ranges over all reactants
    for (int m = 0; m < data.NumReactants; ++m) {
      //double y = Y[ids.Ids[m]];
      //factorForward *= y;
      //if (Ns.Ns[m] > 1)
      //  factorForward *= pow(y, (double)(Ns.Ns[m]-1));
//      if (y > 1.0)
//        printf("WARNING: large Y = %.3E\n", y);
      factorForward *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
    }

    if (alsoDoInverse) {
      // m ranges over all products
      for (int m = data.NumReactants; m < data.NumReactants + data.NumProducts;
          ++m) {
        factorInverse *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
      }
    }

    double contrib = factorForward;
    if (alsoDoInverse) {
      contrib -= factorInverse;

      // if the forward and reverse factors should cancel (but don't due to
      // round-off errors), set contrib to 0
      if (2.0 * fabs(contrib) < 1.0E-8 * (factorForward + factorInverse))
        contrib = 0.0;
    }

    // i ranges over reactants and products
    for (int i = 0; i < data.NumReactants + data.NumProducts; ++i) {
      // in the forward reaction, reactants are destroyed and products are
      // created
      double forwardSign = (i < data.NumReactants ? -1.0 : 1.0);
      double N = (double)Ns.Ns[i];
      (*pYdot)[ids.Ids[i]] += forwardSign * N * contrib;
    }
  }
}

void ReactionLibraryBase::DoAddIndicesOfNonZeroJacobianEntries(
    std::set<std::array<int, 2>> * const indices,
    const bool alsoDoInverse) const {
  for (unsigned int k = 0; k < NumActiveReactions(); ++k) {
    ReactionData_t data = mData[k];
    ReactionIds8_t ids;
    ids.Uint64[0] = mIds0[k].Uint64;
    ids.Uint64[1] = mIds1[k].Uint64;

    for (int j = 0; j < data.NumReactants; ++j) {
      for (int i = 0; i < data.NumProducts; ++i)
        indices->insert( { { ids.Ids[data.NumReactants + i], ids.Ids[j] } });
      for (int i = 0; i < data.NumReactants; ++i)
        indices->insert( { { ids.Ids[i], ids.Ids[j] } });
    }

    if (alsoDoInverse) {
      for (int j = 0; j < data.NumProducts; ++j) {
        for (int i = 0; i < data.NumReactants; ++i)
          indices->insert( { { ids.Ids[i], ids.Ids[data.NumReactants + j] } });
        for (int i = 0; i < data.NumProducts; ++i)
          indices->insert( { { ids.Ids[data.NumReactants + i],
              ids.Ids[data.NumReactants + j] } });
      }
    }
  }
}

void ReactionLibraryBase::DoAddJacobianContributions(
    const std::vector<double>& Y, Jacobian * const pJ,
    const std::vector<double>& rates,
    const std::vector<double>& inverseRates) const {
  bool alsoDoInverse = inverseRates.size() > 0;

  for (unsigned int k = 0; k < NumActiveReactions(); ++k) {
    ReactionData_t data = mData[k];
    ReactionIds8_t ids;
    ids.Uint64[0] = mIds0[k].Uint64;
    ids.Uint64[1] = mIds1[k].Uint64;
    ReactionNs8_t Ns;
    Ns.Uint32[0] = mNs0[k].Uint32;
    Ns.Uint32[1] = mNs1[k].Uint32;

    double rate = rates[k];
    // j ranges over reactants
    for (int j = 0; j < data.NumReactants; ++j) {
      double factor = rate * (double)Ns.Ns[j];
      factor *= pow(Y[ids.Ids[j]], (double)(Ns.Ns[j] - 1));

      // m ranges over all reactants except j
      for (int m = 0; m < data.NumReactants; ++m) {
        if (m == j)
          continue;
        factor *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
      }

      // i ranges over everything
      for (int i = 0; i < data.NumReactants; ++i) {
        (*pJ)(ids.Ids[i], ids.Ids[j]) -= factor * (double)Ns.Ns[i];
      }
      for (int i = data.NumReactants; i < data.NumReactants + data.NumProducts;
          ++i) {
        (*pJ)(ids.Ids[i], ids.Ids[j]) += factor * (double)Ns.Ns[i];
      }
    }

    // inverse reaction
    if (alsoDoInverse) {
      rate = -inverseRates[k];
      // j ranges over reactants
      for (int j = data.NumReactants; j < data.NumReactants + data.NumProducts;
          ++j) {
        double factor = rate * (double)Ns.Ns[j];
        factor *= pow(Y[ids.Ids[j]], (double)(Ns.Ns[j] - 1));

        // m ranges over all reactants except j
        for (int m = data.NumReactants;
            m < data.NumReactants + data.NumProducts; ++m) {
          if (m == j)
            continue;
          factor *= pow(Y[ids.Ids[m]], (double)Ns.Ns[m]);
        }

        // i ranges over everything
        for (int i = 0; i < data.NumReactants; ++i) {
          (*pJ)(ids.Ids[i], ids.Ids[j]) -= factor * (double)Ns.Ns[i];
        }
        for (int i = data.NumReactants;
            i < data.NumReactants + data.NumProducts; ++i) {
          (*pJ)(ids.Ids[i], ids.Ids[j]) += factor * (double)Ns.Ns[i];
        }
      }
    }
  }
}
