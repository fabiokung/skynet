/// \file REACLIBReader.cpp
/// \author jlippuner
/// \since Jul 4, 2014
///
/// \brief
///
///

#include <Utilities/File.hpp>
#include "Reactions/REACLIB.hpp"

#include <array>
#include <cctype>
#include <fstream>


namespace { // unnamed so that these can only be used in this file

struct REACLIBChapter {
  constexpr REACLIBChapter(const int numReactants, const int numProducts) :
      NumReactants(numReactants),
      NumProducts(numProducts) {
  }

  const int NumReactants;
  const int NumProducts;
};

constexpr std::array<REACLIBChapter, 11> REACLIBChapters = { {
    REACLIBChapter(1, 1),
    REACLIBChapter(1, 2),
    REACLIBChapter(1, 3),
    REACLIBChapter(2, 1),
    REACLIBChapter(2, 2),
    REACLIBChapter(2, 3),
    REACLIBChapter(2, 4),
    REACLIBChapter(3, 1),
    REACLIBChapter(3, 2),
    REACLIBChapter(4, 2),
    REACLIBChapter(1, 4)
} };

void MakeReactantNamesAndNs(const std::vector<std::string>& names,
    std::vector<std::string> * const uniqueNames, std::vector<int> * const Ns) {
  for (std::string name : names) {
    // we support the following extension to REACLIB:
    // a nuclide name 'xx#y', where xx (number of characters can vary) are
    // digits and y is either 'n' or 'p', means xx neutron or xx protons
    std::size_t hashPosition = name.find('#');

    // the number of this nuclide in the reaction
    int n = 1;

    if (hashPosition != std::string::npos) {
      // we have a '#', split up the name into number and name
      std::string numberStr = name.substr(0, hashPosition);
      std::string nameStr = name.substr(hashPosition + 1);

      if ((numberStr.size() == 0) || (nameStr.size() == 0))
        throw std::invalid_argument("'" + name + "' is an invalide name");

      n = std::stoi(numberStr);

      if ((nameStr != "n") && (nameStr != "p"))
        throw std::invalid_argument("only neutron and protons can have a "
            "number with a # in front");

      name = nameStr;
    }

    if (n == 0)
      continue;

    // check whether we already have this reactant
    int idx = -1;
    for (unsigned int i = 0; i < uniqueNames->size(); ++i) {
      if (uniqueNames->at(i) == name) {
        idx = i;
        break;
      }
    }

    if (idx >= 0) {
      // we already have this reactant in the list
      Ns->at(idx) += n;
    } else {
      uniqueNames->push_back(name);
      Ns->push_back(n);
    }
  }
}

Reaction REACLIBEntryToReaction(const REACLIBEntry& entry,
    const LeptonMode leptonMode, const NuclideLibrary& nuclib) {
  // fill up reactants vectors, making sure there are no duplicate nuclides
  std::vector<std::string> reactantNames;
  std::vector<int> NsOfReactants;
  MakeReactantNamesAndNs(entry.Reactants(), &reactantNames, &NsOfReactants);

  // fill up products vectors
  std::vector<std::string> productNames;
  std::vector<int> NsOfProducts;
  MakeReactantNamesAndNs(entry.Products(), &productNames, &NsOfProducts);

  bool isDecay = true;

  if (leptonMode == LeptonMode::TreatAllAsDecay) {
    isDecay = true;
  } else if (leptonMode == LeptonMode::TreatAllAsCapture) {
    isDecay = false;
  } else if (leptonMode == LeptonMode::TreatAllAsDecayExceptLabelEC) {
    isDecay = true;
    if (entry.Label() == "ec")
      isDecay = false;
  }

  return Reaction(reactantNames, productNames, NsOfReactants, NsOfProducts,
      entry.IsWeak(), entry.IsInverseRate(), isDecay, entry.Label(), nuclib);
}

} // namespace [unnamed]

REACLIB::REACLIB(const std::string& pathToREACLIBFile) {
  ReadFile(pathToREACLIBFile);
}

bool REACLIB::ReactionContainsUnkownNuclides(const REACLIBEntry& entry,
    const NuclideLibrary& nuclideLibrary) {
  // ignore names that contain a '#'
  for (auto reactant : entry.Reactants()) {
    if (reactant.find('#') != std::string::npos)
      continue;

    if (nuclideLibrary.NuclideIdsVsNames().count(reactant) == 0) {
      //printf("unknown reactant '%s'\n", reactant.c_str());
//      printf("lost %s %s reaction\n", (entry.IsWeak() ? "weak" : "strong"),
//          (entry.IsInverseRate() ? "reverse" : "forward"));
      return true;
    }
  }

  for (auto product : entry.Products()) {
    if (product.find('#') != std::string::npos)
      continue;

    if (nuclideLibrary.NuclideIdsVsNames().count(product) == 0) {
      //printf("unknown product '%s'\n", product.c_str());
//      printf("lost %s %s reaction\n", (entry.IsWeak() ? "weak" : "strong"),
//          (entry.IsInverseRate() ? "reverse" : "forward"));
      return true;
    }
  }

  return false;
}

std::vector<Reaction> REACLIB::GetReactions(const ReactionType type,
    const NuclideLibrary& nuclideLibrary, const LeptonMode leptonMode) const {
  std::vector<Reaction> reactions;

  for (auto entry : mEntries) {
    if (ReactionContainsUnkownNuclides(entry, nuclideLibrary))
      continue;

    if (((type == ReactionType::Weak) && entry.IsWeak())
        || ((type == ReactionType::Strong) && !entry.IsWeak()))
      reactions.push_back(REACLIBEntryToReaction(entry, leptonMode,
          nuclideLibrary));
  }

  return reactions;
}

std::vector<REACLIBCoefficients_t> REACLIB::GetCoefficients(
    const ReactionType type, const NuclideLibrary& nuclideLibrary) const {
  std::vector<REACLIBCoefficients_t> coefficients;

  for (auto entry : mEntries) {
    if (ReactionContainsUnkownNuclides(entry, nuclideLibrary))
      continue;

    if (((type == ReactionType::Weak) && entry.IsWeak())
        || ((type == ReactionType::Strong) && !entry.IsWeak()))
      coefficients.push_back(entry.RateFittingCoefficients());
  }

  return coefficients;
}

std::vector<std::string> REACLIB::GetNuclideNames() const {
  return std::vector<std::string>(mNuclideNames.begin(), mNuclideNames.end());
}

REACLIBEntry REACLIB::ReadEntry(std::istream& is) const {
  // read chapter number
  int chapNum = 0;
  is >> chapNum;

  if ((chapNum <= 0) || (chapNum > 11))
    throw std::invalid_argument("Invalid chapter number: "
        + std::to_string(chapNum));

  REACLIBChapter chap = REACLIBChapters.at(chapNum - 1);

  // got to next line
  std::string line;
  std::getline(is, line);
  if (line.length() != 0)
    throw std::invalid_argument("Expected empty line.");

  // 1st line
  std::getline(is, line);

  // check 5 spaces
  if (line.substr(0, 5) != "     ")
    throw std::invalid_argument("Expected 5 spaces in\n" + line);

  std::size_t pos = 5;

  // read reactants
  std::vector<std::string> reactants;
  for (int i = 0; i < chap.NumReactants; ++i) {
    std::string str = line.substr(pos, 5);
    pos += 5;

    if (str == "     ")
      throw std::invalid_argument("Expected non-empty string in\n" + line);

    // trim leading spaces
    str = str.substr(str.find_first_not_of(" "));
    reactants.push_back(str);
    //printf("%s\n", str.c_str());
  }

  // read products
  std::vector<std::string> products;
  for (int i = 0; i < chap.NumProducts; ++i) {
    std::string str = line.substr(pos, 5);
    pos += 5;

    if (str == "     ")
      throw std::invalid_argument("Expected non-empty string in\n" + line);

    // trim leading spaces
    str = str.substr(str.find_first_not_of(" "));
    products.push_back(str);
    //printf("%s\n", str.c_str());
  }

  // read empty entries
  for (int i = 0; i < 6 - chap.NumReactants - chap.NumProducts; ++i) {
    std::string str = line.substr(pos, 5);
    pos += 5;

    if (str != "     ")
      throw std::invalid_argument("Expected empty string in\n" + line);
  }

  // check 8 spaces
  if (line.substr(pos, 8) != "        ")
    throw std::invalid_argument("Expected 8 spaces in\n" + line);
  pos += 8;

  // read label
  std::string str = line.substr(pos, 4);
  pos += 4;
  str = str.substr(str.find_first_not_of(" "));
  std::string label = str;

  // read rate type
  str = line.substr(pos, 1);
  pos += 1;
  RateType type;
  if ((str == " ") || (str == "n"))
    type = RateType::NonResonant;
  else if (str == "r")
    type = RateType::Resonant;
  else if (str == "w")
    type = RateType::Weak;
  else if (str == "s")
    type = RateType::SpontaneousDecay;
  else
    throw std::invalid_argument("Unrecognized rate type flag '" + str
        + "' in\n" + line);

  // read inverse rate
  str = line.substr(pos, 1);
  pos += 1;
  bool isInverseRate = (str == "v");

  // check 3 spaces
  if (line.substr(pos, 3) != "   ")
    throw std::invalid_argument("Expected 3 spaces in\n" + line);
  pos += 3;

  // read Q value
  str = line.substr(pos, 12);
  pos += 12;
  double q = std::stod(str);

  // check 10 spaces
  if (line.substr(pos, 10) != "          ")
    throw std::invalid_argument("Expected 10 spaces in\n" + line);
  pos += 10;

  // check end of line
  if (pos != line.size())
    throw std::invalid_argument("Expected end of line in\n" + line);

  // read a values
  REACLIBCoefficients_t a;
  for (unsigned int i = 0; i < a.size(); ++i) {
    is >> a[i];
  }

  // finish this line
  std::getline(is, line);
  if (line.length() != 35)
    throw std::invalid_argument("Expected 35 characters in\n" + line);

  return REACLIBEntry(reactants, products, label, type, isInverseRate, q, a);
}

void REACLIB::ReadFile(const std::string& pathToREACLIBFile) {
  if (!File::Exists(pathToREACLIBFile)) {
    throw std::invalid_argument("The REACLIB file '"
        + pathToREACLIBFile + "' does not exist.");
  }

  std::ifstream ifs(pathToREACLIBFile.c_str(), std::ifstream::in);

  while (isalnum(ifs.peek())) {
    auto entry = ReadEntry(ifs);
    mEntries.push_back(entry);
    for (auto& nuclide : entry.Reactants())
      mNuclideNames.insert(nuclide);
    for (auto& nuclide : entry.Products())
      mNuclideNames.insert(nuclide);
  }

  ifs.close();
}
