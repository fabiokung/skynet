/// \file TabulatedNeutrinoSpectrum.cpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <thread>
#include <vector>

#include <gsl/gsl_sf_fermi_dirac.h>

#include "BuildInfo.hpp"
#include "EquationsOfState/EOS.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"
#include "EquationsOfState/NeutrinoDistributionTabulated.hpp"
#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"
#include "EquationsOfState/NeutrinoHistoryTabulated.hpp"
#include "EquationsOfState/TabulatedNeutrinoSpectrum.hpp"
#include "Network/NetworkOptions.hpp"
#include "NuclearData/Nuclide.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/Neutrino.hpp"
#include "Reactions/NeutrinoReactionLibrary.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/FunctionIntegrator.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

namespace {

// the vp-process spectra of Friedland et al. 2026, Table 1
const double TNueMeV = 2.6649665550;
const double EtaNue = 2.1003526050;
const double TNuebarMeV = 3.3917199590;
const double EtaNuebar = 1.4947296660;

// production energy grid: uniform 0.25 MeV to 100 MeV
const double GridStepMeV = 0.25;
const double GridMaxMeV = 100.0;

const double kbMeVPerGK = Constants::BoltzmannConstantInMeVPerGK;

bool failed = false;

void Check(const bool ok, const std::string& what) {
  if (!ok) {
    std::cerr << "FAILED: " << what << std::endl;
    failed = true;
  }
}

void CheckClose(const double got, const double expected, const double tol,
    const std::string& what) {
  const double rel = fabs(got - expected) / fabs(expected);
  if (!(rel < tol)) {
    std::cerr << "FAILED: " << what << ": got " << got << ", expected "
        << expected << ", relative difference " << rel << " >= " << tol
        << std::endl;
    failed = true;
  }
}

double FermiDirac(const double eMeV, const double TMeV, const double eta) {
  return 1.0 / (exp(eMeV / TMeV - eta) + 1.0);
}

std::vector<double> EnergyGrid(const double stepMeV = GridStepMeV,
    const double maxMeV = GridMaxMeV) {
  std::vector<double> energies;
  for (int i = 0; i * stepMeV <= maxMeV + 0.5 * stepMeV; ++i)
    energies.push_back(i * stepMeV);
  return energies;
}

std::vector<double> SampleFermiDirac(const std::vector<double>& energies,
    const double TMeV, const double eta) {
  std::vector<double> values(energies.size());
  for (unsigned int i = 0; i < energies.size(); ++i)
    values[i] = FermiDirac(energies[i], TMeV, eta);
  return values;
}

// int_0^infinity E^n / (exp(E/T - eta) + 1) dE = Gamma(n+1) F_n(eta) T^(n+1)
double AnalyticMoment(const int n, const double TMeV, const double eta) {
  return tgamma(n + 1.0) * gsl_sf_fermi_dirac_int(n, eta) * pow(TMeV, n + 1);
}

// An occupation number nothing like a pinched Fermi-Dirac: a decaying envelope
// with two sharp bumps, so it rises and falls twice (non-monotone). The
// normalization only needs a positive, decaying, integrable shape; this is
// deliberately aggressive to stress the rate integral over a non-FD spectrum.
std::vector<double> NonThermalSpectrum(const std::vector<double>& energies) {
  std::vector<double> f(energies.size());
  for (unsigned int i = 0; i < energies.size(); ++i) {
    const double e = energies[i];
    f[i] = exp(-e / 8.0)
        * (1.0 + 4.0 * exp(-0.5 * pow((e - 12.0) / 3.0, 2))
               + 6.0 * exp(-0.5 * pow((e - 40.0) / 7.0, 2)));
  }
  return f;
}

// BBconst with Gamma(4) F_3(eta) T^4 factored out: converts a third energy
// moment to an energy flux (NeutrinoHistoryTabulated uses the same expression)
double FluxConst() {
  return 0.5 / Constants::Pi * Constants::ErgPerMeV
      / pow(Constants::ReducedPlanckConstantInMeVSec, 3)
      / pow(Constants::SpeedOfLightInCmPerSec, 2);
}

// int_0^inf E^3 (norm f)(E) dE over the handed-out distribution, exactly as the
// reaction library would: to infinity, at the tolerance the distribution asks
// for. No knowledge of the internal norm or moment.
double IntegratedThirdMoment(const std::shared_ptr<NeutrinoDistribution>& dist) {
  FunctionIntegrator integrator(1024, 0.0, dist->IntegrationRelativeError());
  auto f = dist->DistributionFunction(NeutrinoSpecies::NuE);
  return integrator.Integrate([&] (double e) { return e * e * e * f(e); },
      0.0, std::numeric_limits<double>::infinity());
}

// The interpolation error scales with the grid step measured in units of the
// spectrum temperature, so each case gets a grid resolving its own scale. The
// production grid (0.25 MeV to 60 MeV) is the one the vp-process spectra use.
void TestInterpolationAccuracy() {
  struct Case { double T, eta, step, max; };
  const std::vector<Case> cases = {
      { TNueMeV, EtaNue, GridStepMeV, GridMaxMeV },
      { TNuebarMeV, EtaNuebar, GridStepMeV, GridMaxMeV },
      { 1.0, 0.0, 0.1, 100.0 },
      { 10.0, 5.0, 1.0, 300.0 }
  };

  for (auto c : cases) {
    const auto energies = EnergyGrid(c.step, c.max);
    TabulatedNeutrinoSpectrum spec(energies,
        SampleFermiDirac(energies, c.T, c.eta));

    double maxRel = 0.0;
    const double eMax = 2.5 * c.max;
    for (int i = 1; i <= 20000; ++i) {
      const double e = eMax * i / 20000.0;
      maxRel = std::max(maxRel, fabs(spec(e) / FermiDirac(e, c.T, c.eta) - 1.0));
    }

    printf("T = %.4f MeV, eta = %.4f, %.2f MeV grid to %.0f MeV: "
        "max relative interpolation error = %.3e\n", c.T, c.eta, c.step, c.max,
        maxRel);
    Check(maxRel < 1.0e-5, "interpolated Fermi-Dirac spectrum within 1e-5");

    // the moments the normalization and the capture rates depend on
    for (int n : { 2, 3, 5 })
      CheckClose(spec.EnergyMoment(n), AnalyticMoment(n, c.T, c.eta), 1.0e-7,
          "energy moment " + std::to_string(n));
  }
}

// A grid that stops before the spectrum has decayed leaves the log-linear
// extrapolation carrying the flux, which is not a spectrum the table defines.
void TestTooShortGrid() {
  bool threw = false;
  try {
    const auto energies = EnergyGrid(0.5, 120.0);
    TabulatedNeutrinoSpectrum spec(energies,
        SampleFermiDirac(energies, 40.0, 5.0));
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  Check(threw, "energy grid that stops before the spectrum decays is rejected");
}

void TestPositivityAndRinging() {
  const auto energies = EnergyGrid();

  // a hard step is the worst case for a cubic: exponentiating a log-space fit
  // keeps it positive, and the overshoot has to stay bounded by the knots that
  // bracket the query rather than ringing away from them
  std::vector<double> step(energies.size());
  for (unsigned int i = 0; i < energies.size(); ++i)
    step[i] = energies[i] < 20.0 ? 1.0e-2
        : 1.0e-6 * exp(-(energies[i] - 20.0) / 3.0);

  const double maxOvershoot = 10.0;
  bool positive = true;
  bool bounded = true;

  TabulatedNeutrinoSpectrum spec(energies, step);
  for (unsigned int i = 0; i + 1 < energies.size(); ++i) {
    const double lo = std::min(step[i], step[i+1]) / maxOvershoot;
    const double hi = std::max(step[i], step[i+1]) * maxOvershoot;

    for (int k = 1; k < 16; ++k) {
      const double f = spec(energies[i]
          + (energies[i+1] - energies[i]) * k / 16.0);
      if (!(f > 0.0)) positive = false;
      if (f < lo || f > hi) bounded = false;
    }
  }

  Check(positive, "discontinuous spectrum stays positive");
  Check(bounded, "discontinuous spectrum stays within a factor of "
      + std::to_string((int)maxOvershoot) + " of the bracketing knots");

  // a monotone spectrum must stay monotone
  TabulatedNeutrinoSpectrum fd(energies,
      SampleFermiDirac(energies, TNueMeV, EtaNue));
  double previous = fd(0.0);
  bool monotone = true;
  for (int i = 1; i <= 20000; ++i) {
    const double f = fd(GridMaxMeV * i / 20000.0);
    if (f > previous)
      monotone = false;
    previous = f;
  }
  Check(monotone, "monotone spectrum interpolates monotonically");
}

void TestTailExtrapolation() {
  const auto energies = EnergyGrid();
  TabulatedNeutrinoSpectrum spec(energies,
      SampleFermiDirac(energies, TNueMeV, EtaNue));

  for (double e : { 120.0, 200.0, 400.0 })
    CheckClose(spec(e), FermiDirac(e, TNueMeV, EtaNue), 1.0e-4,
        "log-linear tail at " + std::to_string(e) + " MeV");

  // a spectrum that does not decay has no integrable tail
  std::vector<double> growing(energies.size());
  for (unsigned int i = 0; i < energies.size(); ++i)
    growing[i] = 1.0 + energies[i];

  bool threw = false;
  try {
    TabulatedNeutrinoSpectrum bad(energies, growing);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  Check(threw, "non-decaying spectrum is rejected");
}

void TestInvalidInput() {
  const auto energies = EnergyGrid();
  const auto values = SampleFermiDirac(energies, TNueMeV, EtaNue);

  auto rejects = [] (const std::vector<double>& e,
      const std::vector<double>& v, const std::string& what) {
    bool threw = false;
    try {
      TabulatedNeutrinoSpectrum spec(e, v);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    Check(threw, what);
  };

  rejects(energies, std::vector<double>(energies.size() - 1, 1.0),
      "mismatched energy and value counts are rejected");
  rejects({ 0.0, 1.0 }, { 1.0, 0.5 }, "fewer than 3 points are rejected");
  rejects({ 0.0, 200.0, 100.0 }, { 1.0, 0.5, 0.25 },
      "non-monotone energy grid is rejected");

  std::vector<double> withZero = values;
  withZero[10] = 0.0;
  rejects(energies, withZero, "zero spectrum value is rejected");

  std::vector<double> withNegative = values;
  withNegative[10] = -1.0;
  rejects(energies, withNegative, "negative spectrum value is rejected");

  // log f > 700 would overflow exp(); an occupation number never gets there
  std::vector<double> withHuge = values;
  withHuge[10] = 1.0e305;
  rejects(energies, withHuge, "an over-large spectrum value is rejected");

  // steeply-decaying spectrum whose below-grid extrapolation to E = 0 overflows
  // (log f rises by 15 per MeV over the 50 MeV below the grid)
  rejects({ 50.0, 51.0, 52.0 }, { 1.0, exp(-15.0), exp(-30.0) },
      "a below-grid extrapolation that would overflow is rejected");
}

void TestSpeciesHandling() {
  const auto energies = EnergyGrid();
  std::shared_ptr<const TabulatedNeutrinoSpectrum> nue(
      new TabulatedNeutrinoSpectrum(energies,
          SampleFermiDirac(energies, TNueMeV, EtaNue)));

  auto single = NeutrinoDistributionTabulated::Create(
      std::valarray<NeutrinoSpecies>{ NeutrinoSpecies::NuE }, { nue },
      std::valarray<double>{ 1.0 });
  Check(single->DistributionFunction(NeutrinoSpecies::AntiNuE)(10.0) == 0.0,
      "unmatched species has a vanishing distribution function");
  CheckClose(single->DistributionFunction(NeutrinoSpecies::NuE)(10.0),
      FermiDirac(10.0, TNueMeV, EtaNue), 1.0e-5, "single species lookup");

  auto doubled = NeutrinoDistributionTabulated::Create(
      std::valarray<NeutrinoSpecies>{ NeutrinoSpecies::NuE,
          NeutrinoSpecies::NuE }, { nue, nue },
      std::valarray<double>{ 1.0, 2.0 });
  CheckClose(doubled->DistributionFunction(NeutrinoSpecies::NuE)(10.0),
      3.0 * FermiDirac(10.0, TNueMeV, EtaNue), 1.0e-5,
      "repeated species are summed");
}

void TestTailTemperature() {
  const auto energies = EnergyGrid();

  for (auto tEta : { std::make_pair(TNueMeV, EtaNue),
      std::make_pair(TNuebarMeV, EtaNuebar) }) {
    TabulatedNeutrinoSpectrum spec(energies,
        SampleFermiDirac(energies, tEta.first, tEta.second));
    CheckClose(spec.TailTemperatureMeV(), tEta.first, 1.0e-3,
        "tail temperature recovers the Fermi-Dirac temperature");
  }
}

// A tabulated spectrum can't thermalize up to a hotter local fluid the way the
// Fermi-Dirac path does, so it refuses rather than silently diverging.
void TestLocalTemperatureFloor() {
  const auto energies = EnergyGrid();
  std::shared_ptr<const TabulatedNeutrinoSpectrum> spectrum(
      new TabulatedNeutrinoSpectrum(energies,
          SampleFermiDirac(energies, TNueMeV, EtaNue)));
  auto dist = NeutrinoDistributionTabulated::Create(
      std::valarray<NeutrinoSpecies>{ NeutrinoSpecies::NuE }, { spectrum },
      std::valarray<double>{ 1.0 });

  bool threw = false;
  try {
    dist->SetLocalT9(0.5 * TNueMeV / kbMeVPerGK);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  Check(!threw, "a local fluid cooler than the spectrum is accepted");

  threw = false;
  try {
    dist->SetLocalT9(2.0 * TNueMeV / kbMeVPerGK);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  Check(threw, "a local fluid hotter than the spectrum is rejected");
}

// The tabulated history must reproduce the blackbody history, including its
// geometric normalization, when it is handed the same spectrum in tabulated
// form.
void TestHistoryParity() {
  const auto energies = EnergyGrid();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE,
      NeutrinoSpecies::AntiNuE };
  const std::vector<double> T9s = { TNueMeV / kbMeVPerGK,
      TNuebarMeV / kbMeVPerGK };
  const std::vector<double> etas = { EtaNue, EtaNuebar };
  const std::vector<double> lums = { 7.0e51, 5.7e51 };
  const std::vector<std::vector<double>> spectra = {
      SampleFermiDirac(energies, TNueMeV, EtaNue),
      SampleFermiDirac(energies, TNuebarMeV, EtaNuebar) };

  // 1e8 cm puts the neutrinosphere inside the radius, 1e6 cm outside it, so
  // both branches of the extended blackbody normalization are exercised
  for (double radius : { 1.0e8, 1.0e6 }) {
    for (bool pointSource : { true, false }) {
      auto bb = NeutrinoHistoryBlackBody::CreateConstant({ 0.0, 1.0e20 },
          { radius, radius }, species, T9s, etas, lums, false);
      bb.PointSource(pointSource);

      auto tab = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 1.0e20 },
          { radius, radius }, species, energies, spectra, lums, false);
      tab.PointSource(pointSource);

      auto bbDist = bb(1.0);
      auto tabDist = tab(1.0);

      for (auto spec : species) {
        auto bbFunc = bbDist->DistributionFunction(spec);
        auto tabFunc = tabDist->DistributionFunction(spec);

        for (double e : { 0.5, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0 })
          CheckClose(tabFunc(e), bbFunc(e), 1.0e-5,
              "tabulated vs blackbody distribution at " + std::to_string(e)
              + " MeV, r = " + std::to_string(radius) + ", "
              + (pointSource ? "point source" : "sphere"));
      }
    }
  }
}

// Black-box normalization: integrate the distribution the history hands out the
// way the reaction library does -- to infinity, at the distribution's own
// tolerance -- and check it carries the luminosity it was given. Uses an
// arbitrary non-Fermi-Dirac shape, and covers the point source and both
// neutrinosphere branches, at snapshots and interpolated times.
void TestNormalizationByIntegration() {
  // the sharp bumps need a finer grid than the smooth-FD production 0.25 MeV so
  // GSL can certify the distribution's 1e-8 integrating the rates to infinity
  const auto energies = EnergyGrid(0.1, 100.0);
  const auto shape = NonThermalSpectrum(energies);
  const double fluxConst = FluxConst();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE };

  // it really is non-monotone; a Fermi-Dirac occupation never is
  bool rises = false;
  for (unsigned int i = 1; i < shape.size(); ++i)
    if (shape[i] > shape[i-1]) rises = true;
  Check(rises, "the test shape is non-monotone, unlike a Fermi-Dirac");

  // the emitted spectrum's third moment sets the neutrinosphere radius; get it
  // by integrating the raw (unit-norm) shape the same black-box way
  auto rawDist = NeutrinoDistributionTabulated::Create(
      std::valarray<NeutrinoSpecies>{ NeutrinoSpecies::NuE },
      { std::shared_ptr<const TabulatedNeutrinoSpectrum>(
          new TabulatedNeutrinoSpectrum(energies, shape)) },
      std::valarray<double>{ 1.0 });
  const double m3raw = IntegratedThirdMoment(rawDist);

  // the recovered luminosity can only be as accurate as the integration allows
  const double tol = 10.0 * rawDist->IntegrationRelativeError();

  const double L = 7.0e51;
  const double rnu2 = L / (fluxConst * m3raw);
  const double rnu = sqrt(rnu2);

  // point source: L = 4 fluxConst r^2 int E^3 (norm f) dE, at any radius
  for (double radius : { 1.0e8, 3.0e7 }) {
    auto hist = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 1.0e20 },
        { radius, radius }, species, energies, { shape }, { L }, false);
    hist.PointSource(true);
    CheckClose(
        4.0 * fluxConst * IntegratedThirdMoment(hist(1.0)) * radius * radius, L,
        tol, "point-source normalization recovers L (non-FD shape, r = "
        + std::to_string(radius) + ")");
  }

  // neutrinosphere: the handed-out moment is the emitted moment times the
  // geometric dilution W. Measure W as the ratio of the integrated normalized
  // and raw moments and check it against the far-field (r > rnu) and decoupling
  // (r < rnu) formulas.
  for (auto rBranch : { std::make_pair(3.0 * rnu, true),
      std::make_pair(0.5 * rnu, false) }) {
    const double radius = rBranch.first;
    auto hist = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 1.0e20 },
        { radius, radius }, species, energies, { shape }, { L }, false);
    const double wMeasured = IntegratedThirdMoment(hist(1.0)) / m3raw;
    const double wExpected = rBranch.second
        ? 0.5 - 0.5 * sqrt(1.0 - rnu2 / (radius * radius))
        : 1.0 - 0.5 * exp(-6.0 * (rnu2 / (radius * radius) - 1.0));
    CheckClose(wMeasured, wExpected, tol,
        "neutrinosphere dilution from integrated moments (non-FD shape, "
        + std::string(rBranch.second ? "far field" : "decoupling") + ")");
  }

  // time-dependent L(t) and r(t), both interpolated linearly in time and both
  // varying. The check recomputes L(t) and r(t) independently, so it catches an
  // error in either interpolation. Query times hit the snapshots (L and r at
  // knots) and points between them (L and r interpolated); with a shared time
  // axis those are the two reachable knot/interpolated combinations.
  const std::vector<double> times = { 0.0, 1.0, 2.0 };
  const std::vector<double> Ls = { 7.0e51, 5.0e51, 3.0e51 };
  const std::vector<double> radii = { 1.0e8, 1.5e8, 2.0e8 };
  std::vector<std::vector<std::vector<double>>> spectra;
  std::vector<std::vector<double>> lums;
  for (unsigned int t = 0; t < times.size(); ++t) {
    spectra.push_back({ shape });
    lums.push_back({ Ls[t] });
  }
  auto hist = NeutrinoHistoryTabulated::CreateTimeDependent(times, radii,
      species, energies, spectra, lums, false);
  hist.PointSource(true);

  auto lerp = [] (const std::vector<double>& v, double t) {
    return t <= 1.0 ? v[0] + (v[1] - v[0]) * t : v[1] + (v[2] - v[1]) * (t - 1.0);
  };
  for (double t : { 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0 }) {
    const double rt = lerp(radii, t);
    CheckClose(4.0 * fluxConst * IntegratedThirdMoment(hist(t)) * rt * rt,
        lerp(Ls, t), tol, "time-dependent normalization recovers L(t) at t = "
        + std::to_string(t));
  }
}

// Snapshots are interpolated log-linearly in time, so the spectrum halfway
// between two snapshots is their geometric mean.
void TestTimeInterpolation() {
  const auto energies = EnergyGrid();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE };
  const std::vector<std::vector<std::vector<double>>> spectra = {
      { SampleFermiDirac(energies, TNueMeV, EtaNue) },
      { SampleFermiDirac(energies, TNuebarMeV, EtaNuebar) } };
  const std::vector<std::vector<double>> lums = { { 1.0e52 }, { 1.0e52 } };

  auto hist = NeutrinoHistoryTabulated::CreateTimeDependent({ 0.0, 2.0 },
      { 1.0e8, 1.0e8 }, species, energies, spectra, lums, false);
  hist.PointSource(true);

  auto atStart = hist(0.0)->DistributionFunction(NeutrinoSpecies::NuE);
  auto atEnd = hist(2.0)->DistributionFunction(NeutrinoSpecies::NuE);
  auto atMiddle = hist(1.0)->DistributionFunction(NeutrinoSpecies::NuE);

  // the normalization is set by the third moment of the blended shape, so
  // compare shapes relative to a reference energy
  const double eRef = 10.0;
  for (double e : { 1.0, 5.0, 20.0, 40.0 }) {
    const double blended = sqrt(atStart(e) / atStart(eRef)
        * atEnd(e) / atEnd(eRef));
    CheckClose(atMiddle(e) / atMiddle(eRef), blended, 1.0e-12,
        "log-linear blend in time at " + std::to_string(e) + " MeV");
  }

  Check(hist(1.0) == hist(1.0),
      "repeated evaluation at the same time is memoized");
  Check(hist(1.0) != hist(0.5), "a different time rebuilds the distribution");
}

// A non-integrable snapshot must fail when the history is built, not lazily
// during evolution where the reaction library would swallow the error.
void TestEagerValidation() {
  const auto energies = EnergyGrid();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE };
  std::vector<double> flat(energies.size(), 1.0); // does not decay

  bool threw = false;
  try {
    NeutrinoHistoryTabulated::CreateTimeDependent({ 0.0, 1.0 },
        { 1.0e8, 1.0e8 }, species, energies, { { flat }, { flat } },
        { { 1.0e52 }, { 1.0e52 } }, false);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  Check(threw, "a non-decaying snapshot is rejected when the history is built");
}

// The memo behind operator() is the only shared mutable state, so concurrent
// evaluation has to agree with the single-threaded answer.
void TestConcurrentEvaluation() {
  const auto energies = EnergyGrid();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE,
      NeutrinoSpecies::AntiNuE };
  const std::vector<std::vector<double>> spectra = {
      SampleFermiDirac(energies, TNueMeV, EtaNue),
      SampleFermiDirac(energies, TNuebarMeV, EtaNuebar) };

  auto hist = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 10.0 },
      { 1.0e8, 1.0e8 }, species, energies, spectra, { 7.0e51, 5.7e51 }, false);
  hist.PointSource(true);

  std::vector<double> times;
  for (int i = 0; i < 64; ++i)
    times.push_back(10.0 * i / 64.0);

  std::vector<double> expected(times.size());
  for (unsigned int i = 0; i < times.size(); ++i)
    expected[i] = hist(times[i])->DistributionFunction(
        NeutrinoSpecies::NuE)(10.0);

  std::vector<std::thread> threads;
  std::vector<int> mismatches(8, 0);

  for (int t = 0; t < 8; ++t) {
    threads.emplace_back([&, t] () {
      for (int pass = 0; pass < 50; ++pass) {
        for (unsigned int i = 0; i < times.size(); ++i) {
          const double got = hist(times[i])->DistributionFunction(
              NeutrinoSpecies::NuE)(10.0);
          if (got != expected[i])
            ++mismatches[t];
        }
      }
    });
  }

  for (auto& thread : threads)
    thread.join();

  int total = 0;
  for (auto m : mismatches)
    total += m;

  Check(total == 0, "concurrent evaluation agrees with single-threaded "
      "evaluation (" + std::to_string(total) + " mismatches)");
}

// The comparison that matters: the rates the network actually integrates.
void TestRateParity() {
  const double deltaNP = 1.293;
  std::array<double, 24> partitionFunction;
  partitionFunction.fill(1.0);
  std::vector<Nuclide> nuclides = {
      Nuclide(0, 1, 0.0, 0.5, partitionFunction, "n"),
      Nuclide(1, 1, -deltaNP, 0.5, partitionFunction, "p")
  };
  NuclideLibrary nuclib(nuclides, "neutron proton");

  NetworkOptions opts;
  Neutrino lib(SkyNetRoot + "/data/neutrino_reactions.dat", nuclib);

  const auto energies = EnergyGrid();
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE,
      NeutrinoSpecies::AntiNuE };
  const std::vector<double> T9s = { TNueMeV / kbMeVPerGK,
      TNuebarMeV / kbMeVPerGK };
  const std::vector<double> etas = { EtaNue, EtaNuebar };
  const std::vector<double> lums = { 7.0e51, 5.7e51 };
  const std::vector<std::vector<double>> spectra = {
      SampleFermiDirac(energies, TNueMeV, EtaNue),
      SampleFermiDirac(energies, TNuebarMeV, EtaNuebar) };
  const double radius = 1.0e8;

  auto bb = NeutrinoHistoryBlackBody::CreateConstant({ 0.0, 1.0e20 },
      { radius, radius }, species, T9s, etas, lums, false);
  bb.PointSource(true);
  auto tab = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 1.0e20 },
      { radius, radius }, species, energies, spectra, lums, false);
  tab.PointSource(true);

  const double TGK = 5.0;
  const double rhoGCC = 1.0e6;
  std::vector<double> partitionFunctions(nuclib.NumNuclides(), 1.0);

  for (auto mode : { NeutrinoCorrectionMode::None,
      NeutrinoCorrectionMode::WeakMagnetism,
      NeutrinoCorrectionMode::WeakMagnetismExact }) {
    NeutrinoReactionLibrary bbLib(lib, "Blackbody", nuclib, opts, false, true,
        false, mode);
    NeutrinoReactionLibrary tabLib(lib, "Tabulated", nuclib, opts, false, true,
        false, mode);

    ThermodynamicState bbState(TGK, rhoGCC, 0.0, 0.0, 1.0, 0.5, 0.0, bb(1.0));
    ThermodynamicState tabState(TGK, rhoGCC, 0.0, 0.0, 1.0, 0.5, 0.0, tab(1.0));

    bbLib.CalculateRates(bbState, partitionFunctions, opts.RateExpArgumentCap,
        nullptr, nullptr);
    tabLib.CalculateRates(tabState, partitionFunctions, opts.RateExpArgumentCap,
        nullptr, nullptr);

    double maxRel = 0.0;
    for (unsigned int i = 0; i < bbLib.InverseRates().size(); ++i) {
      CheckClose(tabLib.InverseRates()[i], bbLib.InverseRates()[i], 1.0e-7,
          "neutrino capture rate " + std::to_string(i));
      CheckClose(tabLib.Rates()[i], bbLib.Rates()[i], 1.0e-7,
          "electron capture rate " + std::to_string(i));
      maxRel = std::max(maxRel, fabs(tabLib.InverseRates()[i]
          / bbLib.InverseRates()[i] - 1.0));
    }
    printf("correction mode %d: max relative rate difference = %.3e\n",
        (int)mode, maxRel);
  }
}

// A grid too coarse for GSL to integrate the rates must fail loudly, not let
// the reaction library swallow the error and zero the rate.
void TestTooCoarseGridIsFatal() {
  const double deltaNP = 1.293;
  std::array<double, 24> partitionFunction;
  partitionFunction.fill(1.0);
  std::vector<Nuclide> nuclides = {
      Nuclide(0, 1, 0.0, 0.5, partitionFunction, "n"),
      Nuclide(1, 1, -deltaNP, 0.5, partitionFunction, "p")
  };
  NuclideLibrary nuclib(nuclides, "neutron proton");
  NetworkOptions opts;
  Neutrino lib(SkyNetRoot + "/data/neutrino_reactions.dat", nuclib);
  NeutrinoReactionLibrary rl(lib, "Tabulated", nuclib, opts, false, true, false,
      NeutrinoCorrectionMode::None);
  std::vector<double> pf(nuclib.NumNuclides(), 1.0);
  const std::vector<NeutrinoSpecies> species = { NeutrinoSpecies::NuE,
      NeutrinoSpecies::AntiNuE };

  // build the distribution through the history so it carries the point-source
  // normalization the rate library actually integrates (the raw occupation
  // scale would itself change how the integrand converges)
  auto rateIntegrationThrows = [&] (const std::vector<double>& grid,
      const std::vector<double>& nue, const std::vector<double>& nuebar) {
    auto hist = NeutrinoHistoryTabulated::CreateConstant({ 0.0, 1.0e20 },
        { 1.0e8, 1.0e8 }, species, grid, { nue, nuebar }, { 7.0e51, 5.7e51 },
        false);
    hist.PointSource(true);
    ThermodynamicState ts(5.0, 1.0e6, 0.0, 0.0, 1.0, 0.5, 0.0, hist(1.0));
    try {
      rl.CalculateRates(ts, pf, opts.RateExpArgumentCap, nullptr, nullptr);
    } catch (const std::exception&) {
      return true;
    }
    return false;
  };

  // same bumpy shape on both grids, so the difference is only the resolution
  const auto coarse = EnergyGrid(0.5, 100.0);
  Check(rateIntegrationThrows(coarse, NonThermalSpectrum(coarse),
      NonThermalSpectrum(coarse)),
      "a bumpy spectrum on a grid too coarse to integrate the rates fails "
      "loudly instead of silently zeroing the rate");

  const auto fine = EnergyGrid(0.05, 100.0);
  Check(!rateIntegrationThrows(fine, NonThermalSpectrum(fine),
      NonThermalSpectrum(fine)),
      "the same spectrum on a fine grid integrates cleanly");

  const auto production = EnergyGrid();
  Check(!rateIntegrationThrows(production,
      SampleFermiDirac(production, TNueMeV, EtaNue),
      SampleFermiDirac(production, TNuebarMeV, EtaNuebar)),
      "a smooth Fermi-Dirac on the production grid is not rejected");
}

} // namespace

int main(int, char**) {
  FloatingPointExceptions::Enable();

  TestInterpolationAccuracy();
  TestTooShortGrid();
  TestPositivityAndRinging();
  TestTailExtrapolation();
  TestInvalidInput();
  TestSpeciesHandling();
  TestTailTemperature();
  TestLocalTemperatureFloor();
  TestHistoryParity();
  TestNormalizationByIntegration();
  TestTimeInterpolation();
  TestConcurrentEvaluation();
  TestEagerValidation();
  TestRateParity();
  TestTooCoarseGridIsFatal();

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
