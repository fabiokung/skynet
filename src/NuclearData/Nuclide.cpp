/// \file Isotope.cpp
/// \author jlippuner
/// \since Jul 2, 2014
///
/// \brief
///
///

#include <NuclearData/Nuclide.hpp>

#include <cassert>
#include <cctype>
#include <cmath>
#include <stdexcept>
#include <string>

#include "NuclearData/ElementNames.hpp"

constexpr double Nuclide::PartitionFunctionT9[24];
constexpr std::array<double, 24> Nuclide::PartitionFunctionT9_SWIG;

Nuclide::Nuclide(const int z, const int a, const double massExcess,
    const double groundStateSpin,
    const std::array<double, 24>& partitionFunctionLog10,
    const std::string& name) :
    mZ(z),
    mA(a),
    mName(name),
    mMassExcess(massExcess),
    mGroundStateSpin(groundStateSpin),
    mPartitionFunctionLog10(partitionFunctionLog10) {
  assert(mA > 0);
  assert(mZ >= 0);

  // make name
  if (mName.length() > 0)
    return;
  else
    mName = Nuclide::GetName(mZ, mA);
}

std::string Nuclide::GetName(const int Z, const int A) {
  // handle (Z = 0, 1 in a special way)
  if (Z <= 1) {
    if (Z == 0) {
      if (A == 1)
        return "n";
      else
        throw std::invalid_argument("Unexpected A value " + std::to_string(A)
            + " for Z = 0");
    } else if (Z == 1) {
      if (A == 1)
        return "p";
      else if (A == 2)
        return "d";
      else if (A == 3)
        return "t";
      else
        throw std::invalid_argument("Unexpected A value " + std::to_string(A)
            + " for Z = 1");
    } else {
      throw std::invalid_argument("Unexpected Z value " + std::to_string(Z));
    }
  } else {
    try {
      return ElementNameVsZ.at(Z) + std::to_string(A);
    } catch (...) {
      throw std::invalid_argument("Unknown name for Z value "
          + std::to_string(Z));
    }
  }
}

std::array<int, 2> Nuclide::GetAAndZ(const std::string& name) {
  std::string element = "";
  std::string number = "";
  bool readingElement = true;

  for (unsigned int i = 0; i < name.size(); ++i) {
    char c = tolower(name[i]);

    if (isalpha(c)) {
      if (readingElement) {
        element += c;
      } else {
        throw std::runtime_error("Invalid nuclide name '" + name + "'");
      }
    } else if (isdigit(c)) {
      if (readingElement) {
        if (element.size() == 0)
          throw std::runtime_error("Invalid nuclide name '" + name + "'");
        readingElement = false;
      }
      number += c;
    } else {
      throw std::runtime_error("Invalid nuclide name '" + name + "'");
    }
  }

  if (number == "") {
    if (element == "n")
      return {{1, 0}};
    else if (element == "p")
      return {{1, 1}};
    else if (element == "d")
      return {{2, 1}};
    else if (element == "t")
      return {{3, 1}};
    else
      throw std::invalid_argument("'" + name + "' is not a valid nuclide name");
  } else {
    int a = std::stoi(number);
    int z = -1;

    for (auto pair : ElementNameVsZ) {
      if (pair.second == element) {
        z = pair.first;
        break;
      }
    }

    if (z == -1)
      throw std::invalid_argument("'" + element + "' is not a valid element "
          "name");

    return {{a, z}};
  }

  throw std::invalid_argument("Something went wrong when trying to convert'"
      + name + "' to A and Z");
}

Nuclide Nuclide::CreateDummy(const int A, const std::string& name) {
  return Nuclide(A / 2 + 1, A, 0.0, 0.0, std::array<double, 24>(), name);
}


