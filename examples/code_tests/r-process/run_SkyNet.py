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

  # make trajectory
  t0 = 1.0e-3
  tfin = 1.0e9

  nuclib = NuclideLibrary.CreateFromWinv("winvne_v2.0.dat", nuclides)

  opts = NetworkOptions()
  opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
  opts.MassDeviationThreshold = 1.0E-10
  opts.IsSelfHeating = do_heat
  opts.EnableScreening = do_screen
  opts.DisableStdoutOutput = True

  helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

  strongReactionLibrary = REACLIBReactionLibrary("reaclib",
    ReactionType.Strong, do_inv, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Strong reactions", nuclib, opts, True, True)
  weakReactionLibrary = REACLIBReactionLibrary("reaclib",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Weak reactions", nuclib, opts, True, True)
  symmetricFission = REACLIBReactionLibrary("nfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Symmetric neutron induced fission with 0 neutrons emitted",
    nuclib, opts, True, True)
  spontaneousFission = REACLIBReactionLibrary("sfis",
    ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Spontaneous fission", nuclib, opts, True, True)

  reactionLibraries = [strongReactionLibrary, weakReactionLibrary,
      symmetricFission, spontaneousFission]

  screen = SkyNetScreening(nuclib)
  net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)

  #lib = net.GetNuclideLibrary()
  #fout = open("sunet", "w")
  #for n in lib.Names():
    #fout.write("%5s\n" % n)
  #fout.close()

  dat = np.loadtxt("traj_skynet")
  density_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,2], True)
  temperature_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,1], True)

  Ye0 = 0.07
  Temp0 = dat[0,1]
  tfinal = 5.0e8

  # initial abundance
  #nseResult = NSE.CalcFromTemperatureAndDensity(Temp0, density_vs_time(t0),
      #Ye0, net.GetNuclideLibrary())
  #np.savetxt("init_Y", nseResult.Y())

  initY = np.loadtxt("init_Y")

  if (do_heat):
    output = net.EvolveSelfHeatingWithInitialTemperature(initY, t0, tfinal,
        Temp0, density_vs_time, pref)
  else:
    output = net.Evolve(initY, t0, tfinal, temperature_vs_time, density_vs_time, pref)


if __name__ == '__main__':
  num_cores = multiprocessing.cpu_count()
  print "Running with %i worker threads" % num_cores
  pool = multiprocessing.Pool(num_cores)

  args = []

  for inv in [True, False]:
    for scr in [True, False]:
      for heat in [True, False]:
        args.append((inv, scr, heat))

  print(args)

  pool.map_async(run_skynet, args)

  # done submitting jobs
  pool.close()
  pool.join()
