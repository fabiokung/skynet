/// \file NeutrinoDistributionTabulated.hpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONTABULATED_HPP_
#define SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONTABULATED_HPP_

#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/TabulatedNeutrinoSpectrum.hpp"

// Unlike NeutrinoDistributionFermiDirac, this distribution has no temperature
// and so does not apply the local fluid temperature floor set by SetLocalT9();
// the tabulated spectrum is used as given.
class NeutrinoDistributionTabulated : public NeutrinoDistribution {
public:
  static std::shared_ptr<NeutrinoDistribution> Create(
      const std::valarray<NeutrinoSpecies>& species,
      const std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>>&
          spectra,
      const std::valarray<double>& normalization);

  std::function<double(double)> DistributionFunction(
      const NeutrinoSpecies species) const;

  // Unlike the Fermi-Dirac path, a tabulated spectrum cannot thermalize up to
  // a hotter local fluid (the trapping floor in NeutrinoDistributionFermiDirac,
  // T_eff = max(T_nu, T_local)). Rather than silently diverge from that path,
  // throw when the local temperature would exceed the spectrum's own.
  void SetLocalT9(const double T9);

  // The rate integrals run to E -> infinity over this C1 interpolant. GSL
  // certifies 1e-8 as long as the grid resolves the spectrum: 0.25 MeV suffices
  // for a smooth pinched Fermi-Dirac, but a spectrum with sharp structure (an
  // oscillation-mixed nu_e, say) needs a finer grid, or the roundoff at the
  // knots defeats the extrapolation and GSL throws. The analytic path keeps the
  // 1e-12 default.
  double IntegrationRelativeError() const {
    return 1.0e-8;
  }

  // a failed integral means the grid is too coarse to resolve this spectrum;
  // fail loudly rather than let the reaction library zero the rate
  bool IntegrationFailureIsFatal() const {
    return true;
  }

private:
  NeutrinoDistributionTabulated(
      const std::valarray<NeutrinoSpecies>& species,
      const std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>>&
          spectra,
      const std::valarray<double>& normalization);

  std::vector<std::shared_ptr<const TabulatedNeutrinoSpectrum>> mSpectra;
  std::valarray<double> mNormalization;
};

#endif // SRC_EQUATIONSOFSTATE_NEUTRINODISTRIBUTIONTABULATED_HPP_
