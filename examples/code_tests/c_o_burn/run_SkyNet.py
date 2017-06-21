#!/usr/bin/env python

from SkyNet import *
import numpy as np
import sys
import multiprocessing

def run_skynet(args):
  do_inv, do_screen, do_heat = args

  if (do_inv):
    pref = "S_DB_"
  else:
    pref = "S_noDB_"

  if (do_screen):
    pref = pref + "scr_"
  else:
    pref = pref + "noScr_"

  if (do_heat):
    pref = pref + "heat"
  else:
    pref = pref + "noHeat"

  with open("sunet") as f:
    nuclides = [l.strip() for l in f.readlines()]

  nuclib = NuclideLibrary.CreateFromWinv("winvne_v2.0.dat", nuclides)

  opts = NetworkOptions()
  opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
  opts.MassDeviationThreshold = 1.0E-10
  opts.IsSelfHeating = do_heat
  opts.EnableScreening = do_screen
  opts.DisableStdoutOutput = True
  opts.MinDt = 1.0e-22
  opts.NSEEvolutionMinT9 = 20.0

  helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

  strongReactionLibrary = REACLIBReactionLibrary("reaclib",
    ReactionType.Strong, do_inv, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Strong reactions", nuclib, opts, True, True)
  weakReactionLibrary = REACLIBReactionLibrary("reaclib",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Weak reactions", nuclib, opts, True, True)

  reactionLibraries = [strongReactionLibrary, weakReactionLibrary]

  screen = SkyNetScreening(nuclib)
  net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)
  nuclib = net.GetNuclideLibrary()

  fout = open("sunet", "w")
  for n in nuclib.Names():
    fout.write("%5s\n" % n)
  fout.close()

  T = 3
  rho = 1.0e7

  t0 = 0.0
  tfinal = 1.0e2

  Y0 = np.zeros(nuclib.NumNuclides())
  Y0[nuclib.NuclideIdsVsNames()["c12"]] = 0.5 / 12.0
  Y0[nuclib.NuclideIdsVsNames()["o16"]] = 0.5 / 16.0

  temperature_vs_time = ConstantFunction(T)
  density_vs_time = ConstantFunction(rho)

  if (do_heat):
    output = net.EvolveSelfHeatingWithInitialTemperature(Y0, t0, tfinal,
        T, density_vs_time, pref, 1.0E-12)
  else:
    output = net.Evolve(Y0, t0, tfinal, temperature_vs_time,
        density_vs_time, pref, 1.0E-12)


if __name__ == '__main__':
  num_cores = multiprocessing.cpu_count()
  print "Running with %i worker threads" % num_cores
  pool = multiprocessing.Pool(num_cores)

  args = []

  for inv in [True, False]:
    for scr in [True, False]:
#      for heat in [True, False]:
      for heat in [True]:
        args.append((inv, scr, heat))

  print(args)

  pool.map_async(run_skynet, args)

  # done submitting jobs
  pool.close()
  pool.join()
