#!/usr/bin/env python

import sys
import os
import shutil
import multiprocessing
import re

from SkyNet import *
import numpy as np

# location of the history files of the different trajectories
input_dir = "trajectories/"

# location of the final SkyNet output files
output_dir = "output/"

# location of the SkyNet output files while SkyNet is running
working_dir = "work/"

# choose fission reactions
fission_file = "netsu_panov_symmetric_0neut"
#fission_file = "netsu_panov_symmetric_2neut"
#fission_file = "netsu_panov_symmetric_4neut"
#fission_file = "netsu_nfis_Roberts2010rates"

if not os.path.exists(output_dir):
  os.makedirs(output_dir)

if not os.path.exists(working_dir):
  os.makedirs(working_dir)

num_done = multiprocessing.Value('i', 0)


def run_skynet(args):
  global num_done
  outname, total_size = args

  # check if output already exists and exit if so
  if os.path.exists(output_dir + outname + ".h5"):
    num_done.value = num_done.value + 1
    print "[%i/%i] %s: skipped" % (num_done.value, total_size, outname)
    return

  if os.path.exists(output_dir + "FAILED_" + outname + ".h5"):
    num_done.value = num_done.value + 1
    print "[%i/%i] %s: skipped (failed previously)" % (num_done.value, total_size, outname)
    return

  try:
    nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
      + "/data/webnucleo_nuc_v2.0.xml")

    opts = NetworkOptions()
    opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
    opts.MassDeviationThreshold = 1.0E-10
    opts.IsSelfHeating = True
    opts.EnableScreening = True
    opts.DisableStdoutOutput = True

    helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")
    screen = SkyNetScreening(nuclib)

    strongReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType.Strong, True, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Strong reactions", nuclib, opts, True)
    symmetricFission = REACLIBReactionLibrary(SkyNetRoot
      + "/data/netsu_panov_symmetric_0neut", ReactionType.Strong, False,
      LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Symmetric neutron induced fission with 0 neutrons emitted", nuclib, opts,
      False)
    spontaneousFission = REACLIBReactionLibrary(SkyNetRoot +
      "/data/netsu_sfis_Roberts2010rates", ReactionType.Strong, False,
      LeptonMode.TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib, opts,
      False)

    # use only REACLIB weak rates
    weakReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
      ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
      "Weak reactions", nuclib, opts, True)

    net = ReactionNetwork(nuclib, [strongReactionLibrary, weakReactionLibrary,
        symmetricFission, spontaneousFission], helm, screen, opts)

    # read trajectory file
    trajectory_file = input_dir + outname

    # extract times and densities from trajectory_file, also extract initial
    # entropy and electron fraction, for example
    times, densities, entropies = np.loadtxt(trajectory_file, unpack=True, usecols=(0,2,3))
    initial_entropy = entropies[0]

    f = open(trajectory_file, 'r')
    match = re.search("Ye = *([0-9.eE+-]+)", f.read())
    f.close()

    if (match == None):
      num_done.value = num_done.value + 1
      print "[%i/%i] %s: could not read Ye" % (num_done.value, total_size, outname)
      return

    ye = float(match.group(1))

    # if necessary, add time offset
    times = times + 0.003

    density_vs_time = PiecewiseLinearFunction(np.array(times), np.array(densities), False)
    density_profile = PowerLawContinuation(density_vs_time,  -3.0, 2.0)

    net.EvolveSelfHeatingFromNSEWithInitialEntropy(times[0], 1.0E9, initial_entropy,
        density_profile, ye, working_dir + outname)

    # now move the .h5 and .log files from the working directory to the output
    # directory
    shutil.move(working_dir + outname + ".h5", output_dir + outname + ".h5")
    shutil.move(working_dir + outname + ".log", output_dir + outname + ".log")

    num_done.value = num_done.value + 1
    print "[%i/%i] %s: done" % (num_done.value, total_size, outname)
    return

  except Exception as e:
    shutil.move(working_dir + outname + ".h5", output_dir + "FAILED_" + outname + ".h5")
    shutil.move(working_dir + outname + ".log", output_dir + "FAILED_" + outname + ".log")

    num_done.value = num_done.value + 1
    print "[%i/%i] %s: failed with %s " % (num_done.value, total_size, outname, e.message)
    return


if __name__ == '__main__':
  num_cores = multiprocessing.cpu_count()
  print "Running with %i worker threads" % num_cores
  pool = multiprocessing.Pool(num_cores)

  num_done.value = 0

  files = os.listdir(input_dir)

  size = len(files)
  args = [(f, size) for f in files]
  pool.map_async(run_skynet, args)

  # done submitting jobs
  pool.close()
  pool.join()

# monitor progress with
# $ cd work_dir
# $ watch "bash -c 'tail -n 3 *.log'"
