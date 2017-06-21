#!/usr/bin/env python

from SkyNet import *
import numpy as np
import sys
import multiprocessing

def run_skynet(do_nse):
  do_inv = True
  do_screen = False
  do_heat = False

  if (do_nse):
    pref = "S_NSE"
  else:
    pref = "S_noNSE"

  nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml")

  opts = NetworkOptions()
  opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
  opts.MassDeviationThreshold = 1.0E-10
  opts.IsSelfHeating = do_heat
  opts.EnableScreening = do_screen
  opts.DisableStdoutOutput = True
  opts.MinDt = 1.0e-22

  if (do_nse):
    opts.NSEEvolutionMinT9 = 7.0
  else:
    opts.NSEEvolutionMinT9 = 20.0

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
  nuclib = net.GetNuclideLibrary()

  rho = 1.0e8
  ye = 0.1

  t0 = 1.0e-4
  tfinal = 1.0e4

  dat = np.loadtxt("traj")
  temperature_vs_time = PiecewiseLinearFunction(dat[:,0], dat[:,1], True)

  density_vs_time = ConstantFunction(rho)

  nse_opts = NSEOptions()
  if (not do_screen):
    nse_opts.SetDoScreening(False)
    nse = NSE(nuclib, helm, None, nse_opts)
  else:
    nse = NSE(nuclib, helm, screen, nse_opts)

  nseResult = nse.CalcFromTemperatureAndDensity(temperature_vs_time(t0), rho, ye)
  initY = nseResult.Y()

  if (do_heat):
    output = net.EvolveSelfHeatingWithInitialTemperature(initY, t0, tfinal,
        T, density_vs_time, pref, 1.0e-12)
  else:
    output = net.Evolve(initY, t0, tfinal, temperature_vs_time, density_vs_time,
        pref, 1.0e-12)


if __name__ == '__main__':
  num_cores = multiprocessing.cpu_count()
  print "Running with %i worker threads" % num_cores
  pool = multiprocessing.Pool(num_cores)

  args = [True, False]

  pool.map_async(run_skynet, args)

  # done submitting jobs
  pool.close()
  pool.join()
