/// \file InverseRates.cpp
/// \author jlippuner
/// \since Jul 1, 2014
///
/// \brief
///
///

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "BuildInfo.hpp"
#include "Network/NetworkOptions.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/REACLIB.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Utilities/FloatingPointComparison.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

int main(int, char**) {
  FloatingPointExceptions::Enable();

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  REACLIB reaclib(SkyNetRoot + "/data/reaclib");

  auto strongReactions = reaclib.GetReactions(ReactionType::Strong, nuclib,
      LeptonMode::TreatAllAsDecayExceptLabelEC);
  auto strongCoeffs = reaclib.GetCoefficients(ReactionType::Strong, nuclib);
  double rateExpArgumentCap = NetworkOptions().RateExpArgumentCap;

  std::vector<Reaction> reverseStrongReactions;
  std::vector<REACLIBCoefficients_t> reverseStrongCoeffs;
  for (unsigned int i = 0; i < strongReactions.size(); ++i) {
    if (strongReactions[i].IsInverse()) {
      reverseStrongReactions.push_back(strongReactions[i]);
      reverseStrongCoeffs.push_back(strongCoeffs[i]);
    }
  }

  REACLIBReactionLibrary forwardLib(ReactionType::Strong, true,
      "forward strong reactions", "reaclib", strongReactions, strongCoeffs,
      nuclib, NetworkOptions(), false);
  REACLIBReactionLibrary reverseLib(ReactionType::Strong, false,
      "reverse strong reactions", "reaclib", reverseStrongReactions,
      reverseStrongCoeffs, nuclib, NetworkOptions(), false);

  std::multimap<Reaction, int> forwardReactionsIdxMap;
  std::multimap<Reaction, int> reverseReactionsIdxMap;
  std::set<Reaction> forwardReactions;
  std::set<Reaction> reverseReactions;
  std::set<Reaction> weakReactions;

  for (unsigned int i = 0; i < forwardLib.NumAllReactions(); ++i) {
    auto reaction = forwardLib.Reactions()[i];
    forwardReactionsIdxMap.insert({reaction, i});
    forwardReactions.insert(reaction);
  }

  for (unsigned int i = 0; i < reverseLib.NumAllReactions(); ++i) {
    auto reaction = reverseLib.Reactions()[i];
    reverseReactionsIdxMap.insert({reaction, i});
    reverseReactions.insert(reaction);
  }

  printf("Got %lu forward and %lu inverse rates\n", forwardLib.NumAllReactions(),
      reverseLib.NumAllReactions());
  printf("Distinct reactions: %lu forward and %lu inverse\n",
      forwardReactions.size(), reverseReactions.size());

  if (forwardLib.NumAllReactions() != 33338) {
    printf("Wrong number of forward reactions\n");
    return EXIT_FAILURE;
  }

  if (reverseLib.NumAllReactions() != 32995) {
    printf("Wrong number of reverse reactions\n");
    return EXIT_FAILURE;
  }

  if (forwardReactions.size() != 32822) {
    printf("Wrong number of distinct forward reactions\n");
    return EXIT_FAILURE;
  }

  if (reverseReactions.size() != 32479) {
    printf("Wrong number of distinct reverse reactions\n");
    return EXIT_FAILURE;
  }

  unsigned int numInverseReactionsMatched = 0;
  int numBadQ = 0;

  // check Q values
  for (auto& reaction : forwardReactions) {
    double forwardQ = 0.0;
    double reverseQ = 0.0;
    Reaction reverseReaction = reaction.Inverse();

    if (reverseReactionsIdxMap.count(reverseReaction) == 0) {
      //printf("%-38s no reverse rate\n", id.String().c_str());
      continue;
    }

    for (auto iter = forwardReactionsIdxMap.lower_bound(reaction);
        iter != forwardReactionsIdxMap.upper_bound(reaction); ++iter) {
      forwardQ += forwardLib.Qs()[iter->second];
    }

    ++numInverseReactionsMatched;

    for (auto iter = reverseReactionsIdxMap.lower_bound(reverseReaction);
        iter != reverseReactionsIdxMap.upper_bound(reverseReaction); ++iter) {
      reverseQ += reverseLib.Qs()[iter->second];
    }

    double sum = forwardQ + reverseQ;
    if (fcmp(sum, 0.0, 1.0E-10) > 0) {
      printf("bad Q: %-38s %.10e %.10e %.10e\n", reaction.String().c_str(),
          forwardQ, reverseQ, sum);
      ++numBadQ;
    }
  }

  if (numBadQ > 0) {
    printf("Found %i bad Q values.\n", numBadQ);
    return EXIT_FAILURE;
  }

  if (numInverseReactionsMatched != reverseReactions.size()) {
    printf("Could not match all reverse reactions.\n");
    return EXIT_FAILURE;
  }

  // inverse rate
  double T9 = 3.0;
  double rho = 1.0E6;

  ThermodynamicState thermoState(T9, rho, 0.0, 0.0, 0.0, 0.0, 0.0);

  forwardLib.CalculateRates(thermoState,
      nuclib.PartitionFunctionsWithoutSpinTerm(T9), rateExpArgumentCap,
      nullptr, nullptr);
  reverseLib.CalculateRates(thermoState,
      nuclib.PartitionFunctionsWithoutSpinTerm(T9), rateExpArgumentCap,
      nullptr, nullptr);

  NetworkOutput output = NetworkOutput::CreateNew("out.h5", nuclib, false);

  std::vector<std::string> sourcesToConsider{"rath", "ths8", "cf88"};
  int numOk = 0;
  int numBad = 0;

  for (auto& reaction : forwardReactions) {
    Reaction reverseReaction = reaction.Inverse();
    if (reverseReactionsIdxMap.count(reverseReaction) == 0)
      continue;

    int idx = forwardReactionsIdxMap.find(reaction)->second;
    bool skip = true;
    for (auto src : sourcesToConsider) {
      if (src == forwardLib.Labels()[idx]) {
        skip = false;
        break;
      }
    }

    if (skip)
      continue;

    double reverseRateFromForwardReaction = 0.0;
    double reverseRate = 0.0;

    for (auto iter = forwardReactionsIdxMap.lower_bound(reaction);
        iter != forwardReactionsIdxMap.upper_bound(reaction); ++iter) {
      reverseRateFromForwardReaction +=
          forwardLib.InverseRates()[iter->second];
    }

    for (auto iter = reverseReactionsIdxMap.lower_bound(reverseReaction);
        iter != reverseReactionsIdxMap.upper_bound(reverseReaction); ++iter) {
      reverseRate += reverseLib.Rates()[iter->second];
    }

    double fracError = fabs(reverseRateFromForwardReaction - reverseRate)
        / reverseRate;

    if (fracError < 0.1)
      ++numOk;
    else
      ++numBad;
  }

  printf("ok: %i, bad: %i\n", numOk, numBad);

  if (numBad > 11000)
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
