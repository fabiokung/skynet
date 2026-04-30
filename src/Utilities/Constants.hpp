/// \file Constants.hpp
/// \author jlippuner
/// \since Jul 8, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_CONSTANTS_HPP_
#define SKYNET_UTILITIES_CONSTANTS_HPP_

#include <cmath>

// we use the workaround http://stackoverflow.com/a/13406244 to make the
// constants show up under the name Constants in Python instead of the constants
// being global constants

#ifdef SWIG
%rename(Constants) ConstantsStruct;
#endif // SWIG

struct ConstantsStruct {

  // pi
  static constexpr double Pi = 3.1415926535897932385;

  // erg per MeV http://physics.nist.gov/cgi-bin/cuu/Value?tevj
  static constexpr double ErgPerMeV = 1.602176565E-6;

  // c http://physics.nist.gov/cgi-bin/cuu/Value?c
  static constexpr double SpeedOfLightInCmPerSec = 2.99792458E10;

  // hbar http://physics.nist.gov/cgi-bin/cuu/Value?hbarev
  static constexpr double ReducedPlanckConstantInMeVSec = 6.58211928E-22;

  // k_B http://physics.nist.gov/cgi-bin/cuu/Value?tkev
  static constexpr double BoltzmannConstantInMeVPerGK = 8.6173324E-2;

  static constexpr double BoltzmannConstantsInErgPerK =
      BoltzmannConstantInMeVPerGK * 1.0E-9 * ErgPerMeV;

  // N_A http://physics.nist.gov/cgi-bin/cuu/Value?na
  static constexpr double AvogadroConstantInPerGram = 6.02214129E23;

  // 2 pi hbar^2 c^2
  static constexpr double TwoPiHbar2C2 = (2.0 * Pi
      * ReducedPlanckConstantInMeVSec * ReducedPlanckConstantInMeVSec
      * SpeedOfLightInCmPerSec * SpeedOfLightInCmPerSec);

  // 2 pi hbar^2 c^2 N_A^{2/3}
#if defined(__ICC) || defined(__INTEL_COMPILER) || defined(__clang__)
  // pow() is not constexpr in the C++ standard; GCC accepts it as an extension
  // but Intel and Clang do not, so use a precomputed literal instead.
  // Full-precision IEEE 754: pow(6.02214129e23, 2.0/3.0) = 7.13127679759458e15
  static constexpr double TwoPiHbar2C2NA23 = TwoPiHbar2C2 * 7.13127679759458e15;

#else
  static constexpr double TwoPiHbar2C2NA23 = TwoPiHbar2C2
      * pow(AvogadroConstantInPerGram, 2.0 / 3.0);
#endif

  // k_B / (2 pi hbar^2 c^2)
  static constexpr double BoltzmannConstantDivBy2PiHbar2C2InPerMeVPerGKPerCm2 =
      BoltzmannConstantInMeVPerGK / TwoPiHbar2C2;

  // k_B / (2 pi hbar^2 c^2 N_A^{2/3})
  static constexpr double InverseRateFactor = BoltzmannConstantInMeVPerGK
      / TwoPiHbar2C2NA23;

  // m_p http://physics.nist.gov/cgi-bin/cuu/Value?mpc2mev
  static constexpr double ProtonMassInMeV = 938.272046;

  // m_e http://physics.nist.gov/cgi-bin/cuu/Value?mec2mev
  static constexpr double ElectronMassInMeV = 0.510998928;

  // m_p + m_e
  static constexpr double ProtonMassPlusElectronMassInMeV =
      ProtonMassInMeV + ElectronMassInMeV;

  // m_n http://physics.nist.gov/cgi-bin/cuu/Value?mnc2mev
  static constexpr double NeutronMassInMeV = 939.565379;

};

typedef ConstantsStruct Constants;

#endif // SKYNET_UTILITIES_CONSTANTS_HPP_
