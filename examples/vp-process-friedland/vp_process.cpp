/// \file vp_process.cpp
/// \author Fabio Kung
/// \since Apr 2026
///
/// \brief νp-process driver following Friedland et al. 2025/2026.
///
/// Reads a trajectory file (header: "# Ye = X.XX", columns: time[s] T[GK]
/// rho[g/cm3] r[cm]) and evolves the network from NSE with:
///   - REACLIB v2.2 strong + weak reactions
///   - Beard+2017 medium-enhanced triple-alpha (SpecialReactionLibrary)
///   - NeutrinoReactionLibrary from data/neutrino_reactions.dat
///   - Charged-current weak-magnetism/recoil corrections
///   - Exponential L(t) = L0 * exp(-(t - t_ref)/tau_d)
///   - Optional GR blueshift: T_nu *= Phi(r), L_nu *= Phi^4
///
/// Usage:
///   vp_process <trajectory_file> <output_prefix> [options]
///
/// Options:
///   --lbar-ratio R     L_nuebar / L_nue (default 1.0)
///   --L0 L             Initial ν_e luminosity in erg/s (default 7e51)
///   --T-nue T          ν_e temperature in MeV (default 4.0)
///   --T-nuebar T       ν̄_e temperature in MeV (default 5.0)
///   --eta-nue E        ν_e degeneracy parameter (default 0.0)
///   --eta-nuebar E     ν̄_e degeneracy parameter (default 0.0)
///   --tau-d T          Luminosity e-folding time in s (default 3.0)
///   --t-ref T          Reference time for L(t) in s (default 1.0)
///   --M-pns M          PNS mass in solar masses (default 1.4)
///   --no-wm-recoil     Disable weak-magnetism/recoil correction
///   --no-gr            Disable GR blueshift correction
///   --no-alpha         Disable Beard+2017 enhanced triple-alpha
///   --t-end T          End time in s (default 1e9)
///   --rho-slope S      Power-law slope for ρ late-time tail (default -3)
///   --T-slope S        Power-law slope for T late-time tail (default -0.6667)

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

// --- Gravitational constant in CGS and solar mass ---------------------------
static constexpr double GCgs   = 6.67430e-8;   // cm^3 g^-1 s^-2
static constexpr double MsunG  = 1.98892e33;   // g
static constexpr double RRefCm = 5.0e7;        // 500 km reference radius

// Computes Phi(r) = sqrt((1 - 2GM/R_ref) / (1 - 2GM/r))
// Phi > 1 for r < R_ref (blueshift toward PNS), Phi < 1 for r > R_ref.
static double GrBlueshift(const double rCm, const double gmOverC2) {
  const double fRef = 1.0 - 2.0 * gmOverC2 / RRefCm;
  const double fR   = 1.0 - 2.0 * gmOverC2 / rCm;
  if (fR <= 0.0 || fRef <= 0.0) return 1.0;  // inside Schwarzschild radius
  return std::sqrt(fRef / fR);
}

// ---------------------------------------------------------------------------

struct Args {
  std::string TrajFile;
  std::string OutputPrefix;
  double LbarRatio  = 1.0;
  double L0Nue      = 7.0e51;       // erg/s
  double TNueMeV    = 4.0;          // MeV
  double TNuebarMeV = 5.0;          // MeV
  double EtaNue     = 0.0;
  double EtaNuebar  = 0.0;
  double TauD       = 3.0;          // s
  double TRef       = 1.0;          // s, L(t)=L0*exp(-(t-t_ref)/tau_d)
  double MPnsMsun   = 1.4;
  bool   UseWmRecoil = true;
  bool   UseGR      = true;
  bool   UseAlpha   = true;
  double TEnd       = 1.0e9;        // s
  double RhoSlope   = -3.0;         // ρ ∝ t^RhoSlope for late-time tail
  double T9Slope    = -2.0 / 3.0;  // T ∝ t^T9Slope (adiabatic expansion)
  double TailBlend  = 2.0;          // blending window in seconds
};

static Args ParseArgs(int argc, char** argv) {
  if (argc < 3) {
    fprintf(stderr,
        "Usage: %s <trajectory_file> <output_prefix> [options]\n", argv[0]);
    throw std::invalid_argument("too few arguments");
  }

  Args a;
  a.TrajFile     = argv[1];
  a.OutputPrefix = argv[2];

  for (int i = 3; i < argc; ++i) {
    if (!strcmp(argv[i], "--lbar-ratio") && i+1 < argc)
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
    else if (!strcmp(argv[i], "--no-wm-recoil"))
      a.UseWmRecoil = false;
    else if (!strcmp(argv[i], "--no-gr"))
      a.UseGR = false;
    else if (!strcmp(argv[i], "--no-alpha"))
      a.UseAlpha = false;
    else if (!strcmp(argv[i], "--t-end") && i+1 < argc)
      a.TEnd = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--rho-slope") && i+1 < argc)
      a.RhoSlope = std::stod(argv[++i]);
    else if (!strcmp(argv[i], "--T-slope") && i+1 < argc)
      a.T9Slope = std::stod(argv[++i]);
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      throw std::invalid_argument("unknown option");
    }
  }
  return a;
}

// ---------------------------------------------------------------------------
// Trajectory file format:
//   Comment lines start with '#'; one must contain "Ye = <value>".
//   Data lines: time[s]  T[GK]  rho[g/cm3]  r[cm]
// ---------------------------------------------------------------------------
struct Trajectory {
  double Ye = 0.0;
  std::vector<double> Times;
  std::vector<double> TGK;
  std::vector<double> Rho;
  std::vector<double> Radius;
};

static Trajectory ReadTrajectory(const std::string& path) {
  std::ifstream ifs(path);
  if (!ifs.is_open())
    throw std::runtime_error("Cannot open trajectory file: " + path);

  Trajectory traj;
  std::string line;
  bool foundYe = false;

  while (std::getline(ifs, line)) {
    if (line.empty()) continue;

    if (line[0] == '#') {
      // Look for "Ye = X.XX" in header comments
      auto pos = line.find("Ye");
      if (pos != std::string::npos) {
        auto eq = line.find('=', pos);
        if (eq != std::string::npos) {
          traj.Ye = std::stod(line.substr(eq + 1));
          foundYe = true;
        }
      }
      continue;
    }

    std::istringstream iss(line);
    double t, T, rho, r;
    if (!(iss >> t >> T >> rho >> r)) continue;
    traj.Times.push_back(t);
    traj.TGK.push_back(T);
    traj.Rho.push_back(rho);
    traj.Radius.push_back(r);
  }

  if (!foundYe)
    throw std::runtime_error("Trajectory file missing '# Ye = ...' header");
  if (traj.Times.empty())
    throw std::runtime_error("Trajectory file contains no data");

  return traj;
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
  FloatingPointExceptions::Enable();

  Args args = ParseArgs(argc, argv);
  Trajectory traj = ReadTrajectory(args.TrajFile);

  printf("Trajectory: %s  Ye=%.4f  N=%zu  t=[%.3e, %.3e] s\n",
      args.TrajFile.c_str(), traj.Ye,
      traj.Times.size(), traj.Times.front(), traj.Times.back());

  // ---- Network setup -------------------------------------------------------

  auto nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  NetworkOptions opts;
  opts.MaxYChangePerStep                 = 0.1;
  opts.SmallestYUsedForDtCalculation     = 1.0e-4;
  opts.ConvergenceCriterion              = NetworkConvergenceCriterion::Mass;
  opts.SmallestYUsedForErrorCalculation  = 1.0e-20;
  opts.MaxDtChangeMultiplier             = 2.0;
  opts.MinDt                             = 1.0e-16;
  opts.MaxDt                             = 1.0e-3;
  opts.IsSelfHeating                     = false;  // T(t) externally prescribed
  opts.EnableScreening                   = true;

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
      args.UseWmRecoil ? NeutrinoCorrectionMode::WeakMagnetism
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

  // ---- Neutrino history ----------------------------------------------------

  // Convert neutrino temperatures from MeV to GK
  const double kbMeVPerGK   = Constants::BoltzmannConstantInMeVPerGK;
  const double t9NueBase    = args.TNueMeV    / kbMeVPerGK;
  const double t9NuebarBase = args.TNuebarMeV / kbMeVPerGK;

  // GR parameter: GM/c² in cm
  const double c        = Constants::SpeedOfLightInCmPerSec;
  const double gmOverC2 = GCgs * args.MPnsMsun * MsunG / (c * c);

  // Build time-dependent neutrino luminosities and temperatures
  // Two species: NuE (index 0) and AntiNuE (index 1)
  const std::vector<NeutrinoSpecies> nuSpecies = {
      NeutrinoSpecies::NuE, NeutrinoSpecies::AntiNuE
  };

  const std::size_t Nt = traj.Times.size();
  std::vector<double>              radii(Nt);
  std::vector<std::vector<double>> T9s(Nt, std::vector<double>(2));
  std::vector<std::vector<double>> etas(Nt, std::vector<double>(2));
  std::vector<std::vector<double>> lums(Nt, std::vector<double>(2));

  for (std::size_t i = 0; i < Nt; ++i) {
    const double t = traj.Times[i];
    const double r = traj.Radius[i];

    // Luminosity: exponential decay from t_ref
    const double decay = std::exp(-std::max(t - args.TRef, 0.0) / args.TauD);
    double lNue    = args.L0Nue * decay;
    double lNuebar = lNue * args.LbarRatio;

    double t9Nue    = t9NueBase;
    double t9Nuebar = t9NuebarBase;

    if (args.UseGR && r > 0.0) {
      const double phi  = GrBlueshift(r, gmOverC2);
      const double phi4 = phi * phi * phi * phi;
      t9Nue    *= phi;
      t9Nuebar *= phi;
      lNue     *= phi4;
      lNuebar  *= phi4;
    }

    radii[i]   = r;
    T9s[i][0]  = t9Nue;
    T9s[i][1]  = t9Nuebar;
    etas[i][0] = args.EtaNue;
    etas[i][1] = args.EtaNuebar;
    lums[i][0] = lNue;
    lums[i][1] = lNuebar;
  }

  // interpLogSpace=false: etas can be 0, making log-space impossible
  auto nuHist = NeutrinoHistoryBlackBody::CreateTimeDependent(
      traj.Times, radii, nuSpecies, T9s, etas, lums, /*interpLogSpace=*/false);
  nuHist.PointSource(true);

  net.LoadNeutrinoHistory(nuHist.MakeSharedPtr());

  // ---- Trajectory history --------------------------------------------------

  auto hist = TemperatureDensityHistory::CreateFromValues(
      traj.Ye,
      traj.Times.front(),
      traj.Times.back(),
      traj.Times, traj.TGK, traj.Rho);

  // Beyond the trajectory end, extrapolate with power-law tails.
  // PowerLawContinuation blends smoothly from the data into the power law.
  PowerLawContinuation t9WithTail(hist.TemperatureVsTime(),
      args.T9Slope, args.TailBlend);
  PowerLawContinuation rhoWithTail(hist.DensityVsTime(),
      args.RhoSlope, args.TailBlend);

  // ---- Evolve --------------------------------------------------------------

  printf("Evolving from NSE (Ye=%.4f) to t=%.3e s...\n", traj.Ye, args.TEnd);
  printf("  lbar-ratio = %.6f  L0_nue = %.3e erg/s  tau_d = %.1f s\n",
      args.LbarRatio, args.L0Nue, args.TauD);
  printf("  T_nue = %.2f MeV  T_nuebar = %.2f MeV  WM/recoil: %s  GR blueshift: %s\n",
      args.TNueMeV, args.TNuebarMeV, args.UseWmRecoil ? "on" : "off",
      args.UseGR ? "on" : "off");
  printf("  Late-time tails: T ∝ t^%.4f  ρ ∝ t^%.4f\n",
      args.T9Slope, args.RhoSlope);

  net.EvolveFromNSE(
      hist.StartTime(), args.TEnd,
      &t9WithTail, &rhoWithTail,
      hist.Ye(), args.OutputPrefix);

  return EXIT_SUCCESS;
}
