/// \file NuclideLibrary.hpp
/// \author jlippuner
/// \since Jul 16, 2014
///
/// \brief
///
///

#ifndef SKYNET_NUCLEARDATA_NUCLIDELIBRARY_HPP_
#define SKYNET_NUCLEARDATA_NUCLIDELIBRARY_HPP_

#include <map>
#include <string>
#include <unordered_map>
#include <valarray>
#include <vector>

#include "NuclearData/Nuclide.hpp"
#include "Utilities/Interpolators/CubicHermiteInterpolator.hpp"

class NetworkOutput;

class NuclideLibrary {
public:
  /// Create a nuclide library from the given nuclides.
  NuclideLibrary(const std::vector<Nuclide>& nuclides,
      const std::string& source);

  /// Copy constructor, needs to be defined explicitly because of the
  /// std::unique_ptr.
  NuclideLibrary(const NuclideLibrary& other);

  /// Create a nuclide library from a Webnucleo XML file and an optional
  /// explicit list of nuclides that should be in the library. If the optional
  /// list is not given, all nuclides from the Webnucleo XML file are added
  /// to the library.
  static NuclideLibrary CreateFromWebnucleoXML(
      const std::string& webnucleoXMLFile,
      const std::vector<std::string>& explicitListOfNuclideNames =
          std::vector<std::string>());

  /// Create a nuclide library from a Winv file and an optional explicit list of
  /// nuclides that should be in the library. If the optional list is not given,
  /// all nuclides from the Winv file are added to the library.
  static NuclideLibrary CreateFromWinv(const std::string& winvFile,
      const std::vector<std::string>& explicitListOfNuclideNames =
          std::vector<std::string>());

  /// Create a nuclide library from another library by restricting the nuclides
  /// in the original library to the nuclide names given in the restricted list
  /// of nuclide names. Note that no check is made that all the name in the
  /// restricted list actually exist in the original library. Therefore, the
  /// resulting nuclide library will contain no nuclides that are not given in
  /// the list, but it will not contain all nuclides given in the list (it will
  /// contain exactly the intersection of the nuclides in the original library
  /// and the ones whose names are given in the list).
  static NuclideLibrary CreateRestrictedLibrary(
      const NuclideLibrary& originalNuclideLibrary,
      const std::vector<std::string>& restrictedListOfNuclideNames);

  void PrintInfo(NetworkOutput * const pOutput) const;

  const std::vector<double>& PartitionFunctionsWithoutSpinTerm(
      const double T9, const double log10Cap = 100.0) const;

  const std::vector<double>& LogDerivativeOfPartitionFunctionsWithSpinTerm(
      const double T9, const double log10Cap = 100.0) const;

  const std::string& Source() const {
    return mSource;
  }

  int NumNuclides() const {
    return mNuclides.size();
  }

  double NeutronMassInMeV() const {
    return mNeutronMassInMeV;
  }

  double ProtonMassInMeV() const {
    return mProtonMassInMeV;
  }

  const std::vector<int>& Zs() const {
    return mZs;
  }

  const std::vector<int>& As() const {
    return mAs;
  }

  const std::vector<int>& Ns() const {
    return mNs;
  }

  const std::vector<std::string>& Names() const {
    return mNames;
  }

  const std::vector<double>& MassesInMeV() const {
    return mMasses;
  }

  const std::vector<double>& MassExcessesInMeV() const {
    return mMassExcesses;
  }

  const std::vector<double>& BindingEnergiesInMeV() const {
    return mBindingEnergies;
  }

  const std::vector<double>& GroundStateSpinTerms() const {
    return mGroundStateSpinTerms;
  }

  const std::map<std::string, int>& NuclideIdsVsNames() const {
    return mNuclideIdsVsNames;
  }

  const std::vector<Nuclide>& Nuclides() const {
    return mNuclides;
  }

private:
  /// Create a nuclide library from a map that maps nuclides names to Nuclide
  /// class instances. If the explicit list of nuclide names is not empty,
  /// only the nuclides whose names appear in that list are added to the
  /// library, otherwise all nuclides in the map are added to the library.
  NuclideLibrary(const std::unordered_map<std::string, Nuclide>& nuclideMap,
      const std::vector<std::string>& explicitListOfNuclideNames,
      const std::string& source);

  void CalculatePartitionFunctionsAndDerivatives(const double T9,
      const double log10Cap) const;

  void CreatePartitionFunctionInterpolator(const double log10Cap) const;

  std::string mSource;

  // for consistency, if the nuclear data includes data for the neutron and
  // proton, we use that data for the neutron and proton masses, otherwise
  // we use the neutron and proton masses from Utilities/Constants.hpp
  double mNeutronMassInMeV;
  double mProtonMassInMeV;

  // raw nuclear data
  std::vector<int> mZs;
  std::vector<int> mAs;
  std::vector<int> mNs;
  std::vector<std::string> mNames;
  std::vector<double> mMasses; // in MeV
  std::vector<double> mMassExcesses; // in MeV
  std::vector<double> mBindingEnergies; // in MeV
  std::vector<double> mGroundStateSpinTerms;

  // map from the nuclide names to their ids
  std::map<std::string, int> mNuclideIdsVsNames;

  // save the vector of nuclides so that it can be reused to create a restricted
  // nuclide library
  std::vector<Nuclide> mNuclides;

  // computed data, this is mutable because it is used for memoization

  // T9 and log10 cap for which the calculated partition function values and
  // derivatives are bing stored
  mutable double mT9ForWhichPartitionFunctionsWereCalculated;
  mutable double mLog10CapForWhichPartitionFunctionsWereCalcualted;

  // the stored values for the above T9 and log10 cap
  mutable std::vector<double> mPartitionFunctionsWithoutSpinTerm;
  mutable std::vector<double> mLogDerivativeOfPartitionFunctionsWithSpinTerm;

  // the interpolator that interpolates the partition function values and
  // derivatives (each value is a valarray with the partition function for each
  // nuclide), this has to be recreated whenever the log10 cap changes
  mutable std::unique_ptr<CubicHermiteInterpolator<std::valarray<double>>>
      mpPartitionFunctionInterpolator;
};

#endif // SKYNET_NUCLEARDATA_NUCLIDELIBRARY_HPP_
