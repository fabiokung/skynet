/// \file NeutrinoHistoryTabulated.hpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYTABULATED_HPP_
#define SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYTABULATED_HPP_

#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoDistributionTabulated.hpp"

#include <atomic>

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

// Neutrino history driven by a spectrum tabulated on an energy grid, one
// snapshot per time, instead of the (T, eta) pair NeutrinoHistoryBlackBody
// uses. The tabulated spectrum supplies the shape; operator() reuses the
// blackbody flux-dilution factor (point-source 1/r^2, or the solid-angle
// integral of an emitting sphere) verbatim. That factor needs one number from
// the spectrum, the third energy moment, because
//
//   BBconst * T^4 = 1 / (2 pi) * Gamma(4) * F_3(eta) * T^4 * ErgPerMeV / (hbar^3 c^2)
//
// and Gamma(4) * F_3(eta) * T^4 is the third moment of the Fermi-Dirac shape.
// Using the third moment of the tabulated shape reproduces the blackbody
// dilution exactly for a Fermi-Dirac table and generalizes it to an arbitrary
// spectrum. In the point-source case only the shape matters; the sphere
// branches also depend on the moment's absolute scale, so there the spectrum
// must be a true occupation number (see operator()).
class NeutrinoHistoryTabulated : public NeutrinoHistory {
public:
  // spectra are indexed [species][energy], values are unnormalized occupation
  // numbers f(E) sampled at energiesMeV
  static NeutrinoHistoryTabulated CreateConstant(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<double>>& spectra,
      const std::vector<double>& luminosities,
      const bool interpLogSpace) {
    return DoCreateConstant(times, radii, species, energiesMeV, spectra,
        luminosities, interpLogSpace);
  }

  // spectra are indexed [time][species][energy]
  static NeutrinoHistoryTabulated CreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<std::vector<double>>>& spectra,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, radii, species, energiesMeV, spectra,
        luminosities, interpLogSpace);
  }

  // SWIG automatically converts a std::vector of NeutrinoSpecies to a
  // std::vector of ints
  static NeutrinoHistoryTabulated CreateConstantSWIG(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<int>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<double>>& spectra,
      const std::vector<double>& luminosities,
      const bool interpLogSpace) {
    return DoCreateConstant(times, radii,
        ConvertIntVecToNeutrinoSpeciesVec(species), energiesMeV, spectra,
        luminosities, interpLogSpace);
  }

  static NeutrinoHistoryTabulated CreateTimeDependentSWIG(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<int>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<std::vector<double>>>& spectra,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, radii,
        ConvertIntVecToNeutrinoSpeciesVec(species), energiesMeV, spectra,
        luminosities, interpLogSpace);
  }

  std::shared_ptr<NeutrinoDistribution> operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const {
    return
        std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>(
            new NeutrinoHistoryTabulated(*this));
  }

  std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const {
    return std::shared_ptr<NeutrinoHistory>(new NeutrinoHistoryTabulated(*this));
  }

  // Must be set before the history is evaluated: it changes what operator()
  // returns, so it drops the memoized distribution.
  void PointSource(const bool pointSource) {
    mPointSource = pointSource;
    std::atomic_store(&mCache, std::shared_ptr<const Snapshot>());
  }

  void PrintInfo(NetworkOutput * const pOutput) const;

private:
  template<typename T>
  NeutrinoHistoryTabulated(const bool isConstant, const std::vector<T>& species,
      const std::vector<double>& energiesMeV,
      const GeneralPiecewiseLinearFunction<double>& radiusVsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& logSpectrumVsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& lVsTime) :
  NeutrinoHistory(isConstant, species),
  mEnergiesMeV(energiesMeV),
  mRadiusVsTime(radiusVsTime),
  mLogSpectrumVsTime(logSpectrumVsTime),
  mLVsTime(lVsTime), mPointSource(false) {}

  template<typename T>
  static NeutrinoHistoryTabulated DoCreateConstant(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<T>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<double>>& spectra,
      const std::vector<double>& luminosities,
      const bool interpLogSpace);

  template<typename T>
  static NeutrinoHistoryTabulated DoCreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<T>& species,
      const std::vector<double>& energiesMeV,
      const std::vector<std::vector<std::vector<double>>>& spectra,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace);

  std::vector<double> mEnergiesMeV;
  GeneralPiecewiseLinearFunction<double> mRadiusVsTime;
  // log f of each snapshot, [species][energy] flattened species-major. The
  // logs are taken once here rather than on every evaluation, which also means
  // the interpolation between snapshots is log-linear in f and so cannot
  // produce the non-positive values the interpolation in energy rejects.
  GeneralPiecewiseLinearFunction<std::valarray<double>> mLogSpectrumVsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mLVsTime;
  bool mPointSource;

  // The trial step queries the same time ~3x, and building a snapshot re-fits
  // every energy point. Handing out the same distribution twice is safe
  // because the tabulated one never reads the local temperature the consumer
  // writes with SetLocalT9() (it has no temperature to floor). The memo holds
  // the time and distribution in one immutable record published atomically, so
  // a threaded driver can at worst duplicate work, never read a distribution
  // belonging to a different time.
  struct Snapshot {
    double Time;
    std::shared_ptr<NeutrinoDistribution> Distribution;
  };

  mutable std::shared_ptr<const Snapshot> mCache;
};

template<typename T>
NeutrinoHistoryTabulated NeutrinoHistoryTabulated::DoCreateConstant(
    const std::vector<double>& times,
    const std::vector<double>& radii,
    const std::vector<T>& species,
    const std::vector<double>& energiesMeV,
    const std::vector<std::vector<double>>& spectra,
    const std::vector<double>& luminosities,
    const bool interpLogSpace) {
  std::vector<std::vector<std::vector<double>>> spectraVsTime(times.size(),
      spectra);
  std::vector<std::vector<double>> lVsTime(times.size(), luminosities);

  auto hist = DoCreateTimeDependent(times, radii, species, energiesMeV,
      spectraVsTime, lVsTime, interpLogSpace);
  hist.mIsConst = true;

  return hist;
}

template<typename T>
NeutrinoHistoryTabulated NeutrinoHistoryTabulated::DoCreateTimeDependent(
    const std::vector<double>& times,
    const std::vector<double>& radii,
    const std::vector<T>& species,
    const std::vector<double>& energiesMeV,
    const std::vector<std::vector<std::vector<double>>>& spectra,
    const std::vector<std::vector<double>>& luminosities,
    const bool interpLogSpace) {
  if (times.size() != radii.size())
    throw std::invalid_argument("Tabulated neutrino history has "
        + std::to_string(times.size()) + " times but "
        + std::to_string(radii.size()) + " radii");

  if (times.size() != spectra.size())
    throw std::invalid_argument("Tabulated neutrino history has "
        + std::to_string(times.size()) + " times but "
        + std::to_string(spectra.size()) + " spectrum snapshots");

  if (times.size() != luminosities.size())
    throw std::invalid_argument("Tabulated neutrino history has "
        + std::to_string(times.size()) + " times but "
        + std::to_string(luminosities.size()) + " luminosity snapshots");

  const std::size_t numEnergies = energiesMeV.size();

  for (unsigned int i = 0; i < times.size(); ++i) {
    if (species.size() != spectra[i].size())
      throw std::invalid_argument("Tabulated neutrino history snapshot "
          + std::to_string(i) + " has " + std::to_string(spectra[i].size())
          + " spectra but " + std::to_string(species.size()) + " species");

    if (species.size() != luminosities[i].size())
      throw std::invalid_argument("Tabulated neutrino history snapshot "
          + std::to_string(i) + " has "
          + std::to_string(luminosities[i].size()) + " luminosities but "
          + std::to_string(species.size()) + " species");

    for (unsigned int s = 0; s < species.size(); ++s) {
      if (spectra[i][s].size() != numEnergies)
        throw std::invalid_argument("Tabulated neutrino history snapshot "
            + std::to_string(i) + " species " + std::to_string(s) + " has "
            + std::to_string(spectra[i][s].size()) + " values but the energy "
            "grid has " + std::to_string(numEnergies) + " points");
    }
  }

  std::vector<std::valarray<double>> logSpectrumVsTime(times.size(),
      std::valarray<double>(species.size() * numEnergies));
  std::vector<std::valarray<double>> lVsTime(times.size(),
      std::valarray<double>(species.size()));

  for (unsigned int i = 0; i < times.size(); ++i) {
    for (unsigned int s = 0; s < species.size(); ++s) {
      for (unsigned int e = 0; e < numEnergies; ++e) {
        const double f = spectra[i][s][e];
        if (!(f > 0.0) || !std::isfinite(f))
          throw std::invalid_argument("Tabulated neutrino history snapshot "
              + std::to_string(i) + " species " + std::to_string(s)
              + " must be positive and finite, but f("
              + std::to_string(energiesMeV[e]) + " MeV) = "
              + std::to_string(f));

        logSpectrumVsTime[i][s * numEnergies + e] = log(f);
      }

      // Build the snapshot spectrum now so a non-integrable one (does not
      // decay at the top of the grid, or leaves too much flux in the
      // extrapolated tail) fails here at setup rather than lazily inside the
      // network, where NeutrinoReactionLibrary would swallow the integration
      // error and silently zero the rate.
      TabulatedNeutrinoSpectrum(energiesMeV, spectra[i][s]);
    }
    lVsTime[i] = std::valarray<double>(luminosities[i].data(), species.size());
  }

  return NeutrinoHistoryTabulated(false, species, energiesMeV,
      GeneralPiecewiseLinearFunction<double>(times, radii, interpLogSpace),
      // always linear in the stored log f, which is already the log space
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times,
          logSpectrumVsTime, false),
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, lVsTime,
          interpLogSpace));
}

#endif // SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYTABULATED_HPP_
