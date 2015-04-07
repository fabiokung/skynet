/// \file NeutrinoDistributionFermiDirac.hpp
/// \author jlippuner
/// \since Feb 4, 2016
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONFERMIDIRAC_HPP_
#define SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONFERMIDIRAC_HPP_

#include "EquationsOfState/NeutrinoDistribution.hpp"

class NeutrinoDistributionFermiDirac : public NeutrinoDistribution {
public:
  static std::shared_ptr<NeutrinoDistribution> Create(
      const std::valarray<NeutrinoSpecies>& species,
      const std::valarray<double>& T9, const std::valarray<double>& eta);

  static std::shared_ptr<NeutrinoDistribution> Create(
      const std::valarray<NeutrinoSpecies>& species,
      const std::valarray<double>& T9, const std::valarray<double>& eta,
      const std::valarray<double>& normalization);

  std::function<double(double)> DistributionFunction(
      const NeutrinoSpecies species) const;

private:
  NeutrinoDistributionFermiDirac(const std::valarray<NeutrinoSpecies>& species,
      const std::valarray<double>& T9, const std::valarray<double>& eta,
      const std::valarray<double>& normalization);

  std::valarray<double> mTInMeV;  // temperature in MeV
  std::valarray<double> mEta;     // neutrino chemical potential divided by
                                  // temperature (dimensionless)
  std::valarray<double> mNormalization;   // dimensionless normalization factor
};

#endif // SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONFERMIDIRAC_HPP_
