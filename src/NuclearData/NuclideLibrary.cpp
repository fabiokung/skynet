/// \file NuclideLibrary.cpp
/// \author jlippuner
/// \since Jul 16, 2014
///
/// \brief
///
///

#include "NuclearData/NuclideLibrary.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Network/NetworkOutput.hpp"
#include "NuclearData/WebnucleoXmlReader.hpp"
#include "NuclearData/WinvFileReader.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/File.hpp"
#include "Utilities/Exp10.hpp"

namespace { // unnamed so this function can only be used in this file

std::vector<Nuclide> MakeNuclideVector(
    const std::unordered_map<std::string, Nuclide>& nuclideMap,
    const std::vector<std::string>& explicitListOfNuclideNames) {
  std::vector<Nuclide> nuclides;

  if (explicitListOfNuclideNames.size() > 0) {
    for (auto name : explicitListOfNuclideNames) {
      // make sure we have info about this nuclide
      if (nuclideMap.count(name) == 0)
        throw std::invalid_argument("Tried to use unknown nuclide '" + name
            + "'");

      nuclides.push_back(nuclideMap.at(name));
    }
  }
  else {
    for (auto iter : nuclideMap) {
      nuclides.push_back(iter.second);
    }
  }

  return nuclides;
}

bool CompareNuclides(const Nuclide& nucA, const Nuclide& nucB) {
  if (nucA.Z() == nucB.Z()) {
    if (nucA.A() == nucB.A()) {
      return (nucA.Name() < nucB.Name());
    } else {
      return (nucA.A() < nucB.A());
    }
  } else {
    return (nucA.Z() < nucB.Z());
  }
}

} // namespace [unnamed]

NuclideLibrary NuclideLibrary::CreateFromWebnucleoXML(
    const std::string& webnucleoXMLFile,
    const std::vector<std::string>& nuclideNames) {
  return NuclideLibrary(WebnucleoXmlReader::GetAllNuclides(webnucleoXMLFile),
      nuclideNames, webnucleoXMLFile);
}

NuclideLibrary NuclideLibrary::CreateFromWinv(
    const std::string& winvFile,
    const std::vector<std::string>& nuclideNames) {
  return NuclideLibrary(WinvFileReader::GetAllNuclides(winvFile),
      nuclideNames, winvFile);
}

NuclideLibrary NuclideLibrary::CreateRestrictedLibrary(
    const NuclideLibrary& originalNuclideLibrary,
    const std::vector<std::string>& restrictedListOfNuclideNames) {
  std::vector<Nuclide> newNuclides;

  for (unsigned int i = 0; i < restrictedListOfNuclideNames.size(); ++i) {
    std::string name = restrictedListOfNuclideNames[i];
    if (originalNuclideLibrary.NuclideIdsVsNames().count(name) == 0)
      continue;

    int originalId = originalNuclideLibrary.NuclideIdsVsNames().at(name);
    newNuclides.push_back(originalNuclideLibrary.Nuclides()[originalId]);
  }

  return NuclideLibrary(newNuclides, originalNuclideLibrary.Source());
}

//TODO specify order, disallow restriction
NuclideLibrary::NuclideLibrary(const std::vector<Nuclide>& nuclides,
    const std::string& source) :
    mSource(File::MakeAbsolutePath(source)),
    mZs(nuclides.size()),
    mAs(nuclides.size()),
    mNs(nuclides.size()),
    mNames(nuclides.size()),
    mMasses(nuclides.size()),
    mMassExcesses(nuclides.size()),
    mBindingEnergies(nuclides.size()),
    mGroundStateSpinTerms(nuclides.size()),
    mNuclides(nuclides),
    mT9ForWhichPartitionFunctionsWereCalculated(-1.0),
    mLog10CapForWhichPartitionFunctionsWereCalcualted(0.0),
    mPartitionFunctionsWithoutSpinTerm(nuclides.size()),
    mLogDerivativeOfPartitionFunctionsWithSpinTerm(nuclides.size()),
    mpPartitionFunctionInterpolator(nullptr) {
  // sort nuclides
  std::sort(mNuclides.begin(), mNuclides.end(), CompareNuclides);

  // mass = massMultiplier * A + massExcess
  const double massMultiplier = 931.478;

  // for consistency, if the nuclear data includes data for the neutron and
  // proton, we use that data to get the mass excesses of the neutron and
  // proton, otherwise we use the physical constants
  double protonMassExcess = -1.0;
  double neutronMassExcess = -1.0;

  for (unsigned int i = 0; i < mNuclides.size(); ++i) {
    Nuclide nuclide = mNuclides[i];
    mZs[i] = nuclide.Z();
    mAs[i] = nuclide.A();
    mNs[i] = mAs[i] - mZs[i];
    mMassExcesses[i] = nuclide.MassExcess();
    mNames[i] = nuclide.Name();
    mGroundStateSpinTerms[i] = 2.0 * nuclide.GroundStateSpin() + 1.0;

    if ((nuclide.Z() == 0) && (nuclide.A() == 1))
      neutronMassExcess = nuclide.MassExcess();
    else if ((nuclide.Z() == 1) && (nuclide.A() == 1))
      protonMassExcess = nuclide.MassExcess();

    if (mNuclideIdsVsNames.count(nuclide.Name()) != 0)
      throw std::runtime_error("Duplicate nuclide name '" + nuclide.Name()
          + "'");
    mNuclideIdsVsNames[nuclide.Name()] = i;
  }

  if (neutronMassExcess == -1.0) {
    neutronMassExcess = Constants::NeutronMassInMeV - massMultiplier;
  }

  if (protonMassExcess == -1.0) {
    protonMassExcess =
        Constants::ProtonMassPlusElectronMassInMeV - massMultiplier;
  }

  mNeutronMassInMeV = massMultiplier + neutronMassExcess;
  mProtonMassInMeV = massMultiplier + protonMassExcess;

  for (unsigned int i = 0; i < mNuclides.size(); ++i) {
    double massExcess = mMassExcesses[i];
    mBindingEnergies[i] = (double)mZs[i] * protonMassExcess
        + (double)mNs[i] * neutronMassExcess - massExcess;
    mMasses[i] = massMultiplier * (double)mAs[i] + massExcess;
  }
}

NuclideLibrary::NuclideLibrary(const NuclideLibrary& other) :
    mSource(other.mSource),
    mNeutronMassInMeV(other.mNeutronMassInMeV),
    mProtonMassInMeV(other.mProtonMassInMeV),
    mZs(other.mZs),
    mAs(other.mAs),
    mNs(other.mNs),
    mNames(other.mNames),
    mMasses(other.mMasses),
    mMassExcesses(other.mMassExcesses),
    mBindingEnergies(other.mBindingEnergies),
    mGroundStateSpinTerms(other.mGroundStateSpinTerms),
    mNuclideIdsVsNames(other.mNuclideIdsVsNames),
    mNuclides(other.mNuclides),
    mT9ForWhichPartitionFunctionsWereCalculated(
        other.mT9ForWhichPartitionFunctionsWereCalculated),
    mLog10CapForWhichPartitionFunctionsWereCalcualted(
        other.mLog10CapForWhichPartitionFunctionsWereCalcualted),
    mPartitionFunctionsWithoutSpinTerm(
        other.mPartitionFunctionsWithoutSpinTerm),
    mLogDerivativeOfPartitionFunctionsWithSpinTerm(
        other.mLogDerivativeOfPartitionFunctionsWithSpinTerm),
    mpPartitionFunctionInterpolator(nullptr) {
}

NuclideLibrary::NuclideLibrary(
    const std::unordered_map<std::string, Nuclide>& nuclideMap,
    const std::vector<std::string>& explicitListOfNuclideNames,
    const std::string& source) :
    NuclideLibrary(MakeNuclideVector(nuclideMap, explicitListOfNuclideNames),
        source) {}

void NuclideLibrary::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# Nuclide Library\n");
  pOutput->Log("#   Source: %s\n", Source().c_str());
  pOutput->Log("#   Number of nuclides: %i\n", NumNuclides());
  pOutput->Log("#\n");
}

const std::vector<double>& NuclideLibrary::PartitionFunctionsWithoutSpinTerm(
    const double T9, const double log10Cap) const {
  // if we have not already calculated the partition functions for this
  // temperature and log10Cap, calculate it
  if ((T9 != mT9ForWhichPartitionFunctionsWereCalculated)
      || (mLog10CapForWhichPartitionFunctionsWereCalcualted != log10Cap))
    CalculatePartitionFunctionsAndDerivatives(T9, log10Cap);

  return mPartitionFunctionsWithoutSpinTerm;
}

// note that the spin term multiplies the partition function, thus it is an
// additive constant to the log of the partition function, which means it
// goes away once we differentiate with respect to temperature
const std::vector<double>&
NuclideLibrary::LogDerivativeOfPartitionFunctionsWithSpinTerm(
    const double T9, const double log10Cap) const {
  // if we have not already calculated the partition functions for this
  // temperature and log10Cap, calculate it
  if ((T9 != mT9ForWhichPartitionFunctionsWereCalculated)
      || (mLog10CapForWhichPartitionFunctionsWereCalcualted != log10Cap))
    CalculatePartitionFunctionsAndDerivatives(T9, log10Cap);

  return mLogDerivativeOfPartitionFunctionsWithSpinTerm;
}

void NuclideLibrary::CalculatePartitionFunctionsAndDerivatives(const double T9,
    const double log10Cap) const {
  if ((mpPartitionFunctionInterpolator == nullptr)
      || (mLog10CapForWhichPartitionFunctionsWereCalcualted != log10Cap)) {
    CreatePartitionFunctionInterpolator(log10Cap);
    mLog10CapForWhichPartitionFunctionsWereCalcualted = log10Cap;
  }

  mT9ForWhichPartitionFunctionsWereCalculated = T9;

  auto partFunc = mpPartitionFunctionInterpolator->ValueAndFirstDeriv(T9);
  std::valarray<double> values = std::exp(log(10.0) * partFunc.first);
  std::valarray<double> derivs = log(10.0) * T9 * partFunc.second;

  mPartitionFunctionsWithoutSpinTerm.assign(std::begin(values),
      std::end(values));
  mLogDerivativeOfPartitionFunctionsWithSpinTerm.assign(std::begin(derivs),
      std::end(derivs));
}

void NuclideLibrary::CreatePartitionFunctionInterpolator(
    const double log10Cap) const {
  // we need to construct a vector of T9 values, which looks like
  // T9 = [0, 0.99 * minT9, <T9 values from Nuclide::PartitionFunctionT9>,
  //       T9Cap, T9Cap * 1.01, T9Max],
  // where minT9 is the smallest T9 of Nuclide::PartitionFunctionT9, T9Cap is
  // the smallest T9 for which the extrapolated partition function of any
  // nuclide has reached the log10Cap.
  //
  // For each value of T9 we need to construct a valarray with the log10 values
  // of the partition function of each nuclide.
  double T9Max = 10000.0;

  std::size_t numPart = Nuclide::NumPartitionFunctionEntries;
  std::size_t numNuc = mNuclides.size();
  std::vector<double> T9s(numPart + 5);
  std::vector<std::valarray<double>> log10Values(numPart + 5);

  for (unsigned int i = 0; i < numPart; ++i) {
    T9s[i + 2] = Nuclide::PartitionFunctionT9[i];
    log10Values[i + 2].resize(numNuc);
  }

  // loop over all nuclides and do the following:
  // - populate log10Values
  // - find the max T9 for interpolation (above which we extrapolate) for each
  //   nuclide
  // - find the extrapolation coefficients A and B for each nuclide
  // - find the minimum T9 for which the log10 value of all nuclides is below
  //   log10Cap

  std::vector<double> maxT9sForInterpolation(numNuc);

  // we extrapolate as log10 = A + B * T9
  std::vector<double> As(numNuc);
  std::vector<double> Bs(numNuc);

  double minT9Cap = std::numeric_limits<double>::max();

  for (unsigned int n = 0; n < numNuc; ++n) {
    // get partition function log10 values and find the maximum T9 for which
    // the value is below the cap
    unsigned int maxT9Idx = 0;

    for (unsigned int k = 0; k < numPart; ++k) {
      log10Values[k + 2][n] = mNuclides[n].PartitionFunctionLog10()[k];
      if (log10Values[k + 2][n] < Nuclide::PartitionFunctionLog10Cap) {
        maxT9sForInterpolation[n] = Nuclide::PartitionFunctionT9[k];
        maxT9Idx = k;
      }
    }

    // calculate extrapolation coefficients
    if (maxT9Idx == 0)
      throw std::runtime_error("Can't do extrapolation of partition function "
          "for " + mNuclides[n].Name() + " because there are not enough points "
          "below the cap");

    double log10Low = log10Values[maxT9Idx - 1 + 2][n];
    double log10High = log10Values[maxT9Idx + 2][n];
    double T9Low = Nuclide::PartitionFunctionT9[maxT9Idx - 1];
    double T9High = Nuclide::PartitionFunctionT9[maxT9Idx];

    double B = (log10High - log10Low) / (T9High - T9Low);
    double A = log10Low - B * T9Low;

    As[n] = A;
    Bs[n] = B;

    if (B != 0.0) {
      double T9Cap = (log10Cap - A) / B;
      minT9Cap = std::min(minT9Cap, T9Cap);
    }
  }

  if (minT9Cap >= T9Max)
    minT9Cap = 0.98 * T9Max;

  // loop over nuclides again and perform extrapolation where necessary
  T9s[numPart + 2] = minT9Cap;
  log10Values[numPart + 2].resize(numNuc);

  for (unsigned int n = 0; n < numNuc; ++n) {
    double A = As[n];
    double B = Bs[n];
    double maxT9 = maxT9sForInterpolation[n];

    for (unsigned int k = 0; k < numPart; ++k) {
      if (Nuclide::PartitionFunctionT9[k] > maxT9)
        log10Values[k + 2][n] = A + B * Nuclide::PartitionFunctionT9[k];
    }

    log10Values[numPart + 2][n] = A + B * minT9Cap;
  }

  // we extend the partition function as a constant down to T9 = 0
  T9s[0] = 0.0;
  log10Values[0] = log10Values[2];

  T9s[1] = 0.99 * Nuclide::PartitionFunctionT9[0];
  log10Values[1] = log10Values[2];

  // and as a constant up to T9 = T9Max
  T9s[numPart + 3] = minT9Cap * 1.01;
  log10Values[numPart + 3] = log10Values[numPart + 2];

  T9s[numPart + 4] = T9Max;
  log10Values[numPart + 4] = log10Values[numPart + 2];

  mpPartitionFunctionInterpolator =
      std::unique_ptr<CubicHermiteInterpolator<std::valarray<double>>>(
          new CubicHermiteInterpolator<std::valarray<double>>(T9s,
              log10Values));
}
