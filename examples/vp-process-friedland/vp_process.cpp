/// \file vp_process.cpp
/// \author Fabio Kung
/// \since Apr 2026
///
/// \brief νp-process driver following Friedland et al. 2025/2026.
///
/// Usage:
///   vp_process <trajectory_file> <output_prefix> [options]
///
/// Options:
///   --ye Y             Initial electron fraction (default 0.6)
///   --lbar-ratio R     L_nuebar / L_nue (default 1.0)
///   --L0 L             Initial ν_e luminosity in erg/s (default 7e51)
///   --T-nue T          ν_e temperature in MeV (default 2.67)
///   --T-nuebar T       ν̄_e temperature in MeV (default 3.39)
///   --eta-nue E        ν_e degeneracy parameter (default 2.1)
///   --eta-nuebar E     ν̄_e degeneracy parameter (default 1.5)
///   --tau-d T          Luminosity e-folding time in s (default 3.0)
///   --t-ref T          Reference time for L(t) in s (default 1.0)
///   --M-pns M          PNS mass in solar masses (default 1.4)
///   --wm-mode M        Weak-magnetism/recoil correction: none|first|exact
///                      (first = Horowitz Eq.24/25, exact = Eq.22; default first)
///   --no-wm-recoil     Alias for --wm-mode none
///   --no-gr            Disable GR blueshift correction
///   --no-alpha         Disable Beard+2017 enhanced triple-alpha
///   --t-end T          End time in s (default 1e9)
///   --max-dt T         Maximum network time step in s (default 1e8)
///   --rho-slope S      Power-law slope for ρ late-time tail (default -2)
///   --T-slope S        Power-law slope for T late-time tail (default -0.6667)
///   --tail-blend T     Blend span in s into the power-law tail, must be < trajectory span (default: 1/3 of span)
///   --smallest-y-dt Y  Species with Y below this don't constrain the time step (default 1e-6)
///   --max-ychange Y    Max |dY/Y| per step on the dt-controlling species (default 0.1)
///   --newton-crit M    Newton convergence criterion: mass|dyby|both (default mass)

#include "BuildInfo.hpp"

#include <cstdio>
#include <cmath>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "DensityProfiles/PowerLawContinuation.hpp"
#include "EquationsOfState/HelmholtzEOS.hpp"
#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"
#include "EquationsOfState/SkyNetScreening.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/ReactionNetwork.hpp"
#include "Network/TemperatureDensityHistory.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Reactions/NeutrinoReactionLibrary.hpp"
#include "Reactions/REACLIBReactionLibrary.hpp"
#include "Reactions/SpecialReactionLibrary.hpp"
#include "Utilities/Constants.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

#include "AlphaBurning.hpp"

static constexpr double GCgs   = 6.67430e-8;  // cm^3 g^-1 s^-2
static constexpr double MsunG  = 1.98892e33;  // g
static constexpr double RRefCm = 5.0e7;       // 500 km reference radius

// Tracer (SI) -> SkyNet (CGS/GK) unit conversions.
static constexpr double MToCm      = 100.0;   // m -> cm
static constexpr double KgM3ToGCm3 = 1.0e-3;  // kg/m^3 -> g/cm^3
static constexpr double JToMeV     = 1.0e7 / Constants::ErgPerMeV;

// Computes Phi(r) = sqrt((1 - 2GM/R_ref) / (1 - 2GM/r))
// Phi > 1 for r < R_ref (blueshift toward PNS), Phi < 1 for r > R_ref.
static double GrBlueshift(const double rCm, const double gmOverC2) {
  const double fRef = 1.0 - 2.0 * gmOverC2 / RRefCm;
  const double fR   = 1.0 - 2.0 * gmOverC2 / rCm;
  if (fR <= 0.0 || fRef <= 0.0) return 1.0;  // inside Schwarzschild radius
  return std::sqrt(fRef / fR);
}

struct Args {
  std::string TrajFile;
  std::string OutputPrefix;
  double Ye = 0.6;
  double LbarRatio = 1.0;
  double L0Nue = 7.0e51;
  double TNueMeV = 2.67;
  double TNuebarMeV = 3.39;
  double EtaNue = 2.1;
  double EtaNuebar = 1.5;
  double TauD = 3.0;
  double TRef = 1.0; // L(t)=L0*exp(-(t-t_ref)/tau_d)
  double MPnsMsun = 1.4;
  std::string WmMode = "first"; // weak-magnetism correction: none|first|exact
  bool UseGR = true;
  bool UseAlpha = true;
  double TEnd = 1.0e9;
  double MaxDt = 1.0e8;
  double SmallestYForDt = 1.0e-6; // species below this don't constrain dt; SkyNet's default
  double MaxYChange = 0.1; // max |dY/Y| per step on dt-controlling species
  std::string NewtonCrit = "mass"; // Newton convergence: mass|dyby|both
  double RhoSlope = -2.0; // ρ ∝ t^RhoSlope for late-time tail
  double T9Slope = -2.0 / 3.0; // T ∝ t^T9Slope (adiabatic expansion)
  double TailBlend = 0.0; // 0 => auto: 1/3 of the trajectory time span
};

static Args ParseArgs(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr,
        "Usage: %s <trajectory_file> <output_prefix> [options]\n", argv[0]);
    throw std::invalid_argument("too few arguments");
  }

  Args a;
  a.TrajFile = argv[1];
  a.OutputPrefix = argv[2];

  for (int i = 3; i < argc; ++i) {
    if (!strcmp(argv[i], "--ye") && i+1 < argc)
      a.Ye = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--lbar-ratio") && i+1 < argc)
      a.LbarRatio = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--L0") && i+1 < argc)
      a.L0Nue = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--T-nue") && i+1 < argc)
      a.TNueMeV = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--T-nuebar") && i+1 < argc)
      a.TNuebarMeV = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--eta-nue") && i+1 < argc)
      a.EtaNue = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--eta-nuebar") && i+1 < argc)
      a.EtaNuebar = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--tau-d") && i+1 < argc)
      a.TauD = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--t-ref") && i+1 < argc)
      a.TRef = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--M-pns") && i+1 < argc)
      a.MPnsMsun = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--wm-mode") && i+1 < argc) {
      a.WmMode = argv[++i];
      if (a.WmMode != "none" && a.WmMode != "first" && a.WmMode != "exact")
        throw std::invalid_argument("--wm-mode must be none|first|exact");
    }
    else if (!strcmp(argv[i], "--no-wm-recoil"))
      a.WmMode = "none";
    else if (!strcmp(argv[i], "--no-gr"))
      a.UseGR = false;
    else if (!strcmp(argv[i], "--no-alpha"))
      a.UseAlpha = false;
    else if (!strcmp(argv[i], "--t-end") && i+1 < argc)
      a.TEnd = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--max-dt") && i+1 < argc)
      a.MaxDt = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--rho-slope") && i+1 < argc)
      a.RhoSlope = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--T-slope") && i+1 < argc)
      a.T9Slope = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--tail-blend") && i+1 < argc)
      a.TailBlend = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--smallest-y-dt") && i+1 < argc)
      a.SmallestYForDt = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--max-ychange") && i+1 < argc)
      a.MaxYChange = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--newton-crit") && i+1 < argc) {
      a.NewtonCrit = argv[++i];
      if (a.NewtonCrit != "mass" && a.NewtonCrit != "dyby" &&
          a.NewtonCrit != "both")
        throw std::invalid_argument("--newton-crit must be mass|dyby|both");
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      throw std::invalid_argument("unknown option");
    }
  }
  return a;
}

// Tracer-particle trajectory in SI units, one row per time, whitespace-separated.
// Lines starting with '#' are comments. Six columns, in order:
//   time[s]  radius[m]  velocity[m/s]  kT[J]  entropy[k_B/baryon]  rho[kg/m^3]
// Entropy is read but unused; velocity is kept to coast the radius past the
// trajectory end for the neutrino flux. Ye is not in the file; it is supplied
// via --ye. Stored fields are converted to SkyNet units (GK, cm, cm/s, g/cm^3)
// on read so the rest of the driver works in CGS/GK throughout.
struct Trajectory {
  double Ye = 0.0;
  std::vector<double> Times;
  std::vector<double> TGK;
  std::vector<double> Rho;
  std::vector<double> Radius;
  std::vector<double> Vel;
};

static Trajectory ReadTrajectory(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs.is_open())
    throw std::runtime_error("Cannot open trajectory file: " + path);

  Trajectory traj;
  std::string line;

  while (std::getline(ifs, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;

    std::istringstream iss(line);
    double t, rMeters, vel, kTJoule, entropy, rhoSI;
    if (!(iss >> t >> rMeters >> vel >> kTJoule >> entropy >> rhoSI))
      throw std::runtime_error("Malformed trajectory row in " + path + ": " + line);

    traj.Times.push_back(t);
    traj.TGK.push_back(kTJoule * JToMeV / Constants::BoltzmannConstantInMeVPerGK);
    traj.Rho.push_back(rhoSI * KgM3ToGCm3);
    traj.Radius.push_back(rMeters * MToCm);
    traj.Vel.push_back(vel * MToCm);
  }

  if (traj.Times.empty())
    throw std::runtime_error("Trajectory file contains no data");

  return traj;
}

int main(int argc, char** argv) {
  FloatingPointExceptions::Enable();

  Args args = ParseArgs(argc, argv);
  Trajectory traj = ReadTrajectory(args.TrajFile);
  traj.Ye = args.Ye;

  printf("Trajectory: %s  Ye=%.4f  N=%zu  t=[%.3e, %.3e] s\n",
      args.TrajFile.c_str(), traj.Ye,
      traj.Times.size(), traj.Times.front(), traj.Times.back());

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.MaxYChangePerStep = args.MaxYChange;
  opts.SmallestYUsedForDtCalculation = args.SmallestYForDt;
  opts.ConvergenceCriterion =
      args.NewtonCrit == "dyby" ? NetworkConvergenceCriterion::DeltaYByY
      : args.NewtonCrit == "both" ? NetworkConvergenceCriterion::BothDeltaYByYAndMass
      : NetworkConvergenceCriterion::Mass;
  opts.SmallestYUsedForErrorCalculation = 1.0e-20;
  opts.MaxDtChangeMultiplier = 2.0;
  opts.MinDt = 1.0e-16;
  opts.MaxDt = args.MaxDt;
  opts.IsSelfHeating = false; // T(t) externally prescribed
  opts.EnableScreening = true;

  REACLIBReactionLibrary weakLib(
      SkyNetRoot + "/data/reaclib",
      ReactionType::Weak, false,
      LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, true);

  REACLIBReactionLibrary strongLib(
      SkyNetRoot + "/data/reaclib",
      ReactionType::Strong, true,
      LeptonMode::TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, true);

  // Neutrino reactions (ν_e n ↔ e⁻ p and ν̄_e p ↔ e⁺ n)
  // onlyNuCap=false: include both ν capture and e capture
  // nuHeating=false: T(t) is externally prescribed, not self-consistent
  // includeBeta=true: include β decays
  NeutrinoReactionLibrary nuLib(
      SkyNetRoot + "/data/neutrino_reactions.dat",
      "Neutrino reactions", nuclib, opts,
      /*onlyNuCap=*/false, /*nuHeating=*/false, /*includeBeta=*/true,
      args.WmMode == "exact" ? NeutrinoCorrectionMode::WeakMagnetismExact
          : args.WmMode == "first" ? NeutrinoCorrectionMode::WeakMagnetism
          : NeutrinoCorrectionMode::None);

  // Alpha burning (Beard+2017 medium-enhanced triple-alpha)
  // Be9Fac=1 (nominal), enhancement factors from RunSetup.cpp defaults
  std::unique_ptr<SpecialReactionLibrary> alphaLib;
  if (args.UseAlpha) {
    alphaLib = std::unique_ptr<SpecialReactionLibrary>(
        new SpecialReactionLibrary(
            AlphaBurning::GetRates(nuclib, /*enhanced=*/true,
                /*Be9Fac=*/1.0, /*enhancePEnhance=*/1.0, /*enhanceNEnhance=*/1.0),
            ReactionType::Strong, "Beard+2017 enhanced triple-alpha",
            "Beard+2017", nuclib, opts));
    // Remove alpha-burning reactions from strongLib to avoid double-counting
    strongLib.RemoveReactions(alphaLib->Reactions().AllData());
    printf("Alpha burning: enabled (Beard+2017)\n");
  } else {
    printf("Alpha burning: disabled\n");
  }

  // Assemble reaction library list
  std::vector<const ReactionLibraryBase*> reactionLibs = {
      &weakLib, &strongLib, &nuLib
  };
  if (alphaLib)
    reactionLibs.push_back(alphaLib.get());

  HelmholtzEOS helmEOS(SkyNetRoot + "/data/helm_table.dat");
  SkyNetScreening screen(nuclib);

  ReactionNetwork net(nuclib, reactionLibs, &helmEOS, &screen, opts);

  // neutrino history

  const double kbMeVPerGK = Constants::BoltzmannConstantInMeVPerGK;
  const double t9NueBase = args.TNueMeV / kbMeVPerGK;
  const double t9NuebarBase = args.TNuebarMeV / kbMeVPerGK;

  const double c = Constants::SpeedOfLightInCmPerSec;
  const double gmOverC2 = GCgs * args.MPnsMsun * MsunG / (c * c);

  const std::vector<NeutrinoSpecies> nuSpecies = {
      NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE
  };

  // The neutrino history must span the whole evolution [t0, TEnd], but the tracer
  // ends at traj.Times.back(). L(t), T9 and eta are analytic in t and extend on
  // their own; only the radius needs continuing, because the point-source flux
  // seen by the parcel scales as 1/r^2. By the tracer's end the parcel has left
  // the acceleration region and coasts at its asymptotic wind velocity, so we
  // continue it ballistically, r(t) = r_end + v_end*(t - t_end) -- the outflow
  // behaviour assumed in Friedland et al. (2026). The residual ν flux this far
  // out is a small correction; alternative continuations (holding r fixed, or
  // homologous r ∝ t) barely move the final yields.
  std::vector<double> nuTimes = traj.Times;
  std::vector<double> nuRadii = traj.Radius;
  const double tLast = traj.Times.back();
  const double tGridEnd = args.TEnd * 1.01; // cover steps that probe just past TEnd
  if (tGridEnd > tLast) {
    const double rLast = traj.Radius.back();
    const double vLast = traj.Vel.back();
    const int nExt = 300;
    const double ratio = std::pow(tGridEnd / tLast, 1.0 / nExt);
    double t = tLast;
    for (int k = 0; k < nExt; ++k) {
      t *= ratio;
      nuTimes.push_back(t);
      nuRadii.push_back(rLast + vLast * (t - tLast));
    }
  }

  const std::size_t Nnu = nuTimes.size();
  std::vector<std::vector<double>> T9s(Nnu, std::vector<double>(2));
  std::vector<std::vector<double>> etas(Nnu, std::vector<double>(2));
  std::vector<std::vector<double>> lums(Nnu, std::vector<double>(2));

  for (std::size_t i = 0; i < Nnu; ++i) {
    const double t = nuTimes[i];
    const double r = nuRadii[i];

    const double decay = std::exp(-std::max(t - args.TRef, 0.0) / args.TauD);
    double lNue = args.L0Nue * decay;
    double lNuebar = lNue * args.LbarRatio;
    double t9Nue = t9NueBase;
    double t9Nuebar = t9NuebarBase;

    if (args.UseGR && r > 0.0) {
      const double phi = GrBlueshift(r, gmOverC2);
      const double phi4 = phi * phi * phi * phi;
      t9Nue *= phi;
      t9Nuebar *= phi;
      lNue *= phi4;
      lNuebar *= phi4;
    }

    T9s[i][0] = t9Nue;
    T9s[i][1] = t9Nuebar;
    etas[i][0] = args.EtaNue;
    etas[i][1] = args.EtaNuebar;
    lums[i][0] = lNue;
    lums[i][1] = lNuebar;
  }

  // interpLogSpace=false: etas can be 0, making log-space impossible
  auto nuHist = NeutrinoHistoryBlackBody::CreateTimeDependent(
      nuTimes, nuRadii, nuSpecies, T9s, etas, lums, /*interpLogSpace=*/false);
  nuHist.PointSource(true);

  net.LoadNeutrinoHistory(nuHist.MakeSharedPtr());

  auto hist = TemperatureDensityHistory::CreateFromValues(
      traj.Ye,
      traj.Times.front(),
      traj.Times.back(),
      traj.Times, traj.TGK, traj.Rho);

  // Beyond the trajectory end, extrapolate with power-law tails.
  // PowerLawContinuation blends smoothly from the data into the power law over
  // the last tailBlend seconds; the blend span must be shorter than the data.
  const double trajSpan = traj.Times.back() - traj.Times.front();
  const double tailBlend =
      args.TailBlend > 0.0 ? args.TailBlend : trajSpan / 3.0;
  if (tailBlend >= trajSpan)
    throw std::runtime_error(
        "--tail-blend must be smaller than the trajectory span ("
        + std::to_string(trajSpan) + " s)");

  PowerLawContinuation t9WithTail(hist.TemperatureVsTime(),
      args.T9Slope, tailBlend);
  PowerLawContinuation rhoWithTail(hist.DensityVsTime(),
      args.RhoSlope, tailBlend);

  printf("Evolving from NSE (Ye=%.4f) to t=%.3e s...\n", traj.Ye, args.TEnd);
  printf("  lbar-ratio = %.6f  L0_nue = %.3e erg/s  tau_d = %.1f s\n",
      args.LbarRatio, args.L0Nue, args.TauD);
  printf("  T_nue = %.2f MeV  T_nuebar = %.2f MeV  WM mode: %s  GR blueshift: %s\n",
      args.TNueMeV, args.TNuebarMeV, args.WmMode.c_str(),
      args.UseGR ? "on" : "off");
  printf("  Late-time tails: T ∝ t^%.4f  ρ ∝ t^%.4f\n",
      args.T9Slope, args.RhoSlope);
  printf("  Max dt = %.3e s\n", args.MaxDt);

  net.EvolveFromNSE(
      hist.StartTime(), args.TEnd,
      &t9WithTail, &rhoWithTail,
      hist.Ye(), args.OutputPrefix);

  return EXIT_SUCCESS;
}
