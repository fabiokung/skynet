#!/usr/bin/env python

from SkyNet import *
import numpy as np
import sys

nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

opts = NetworkOptions()
opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
opts.MassDeviationThreshold = 1.0E-10
opts.IsSelfHeating = False
opts.EnableScreening = True
opts.NSEEvolutionMinT9 = 10.0
opts.DisableStdoutOutput = True

helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

strongReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
  ReactionType.Strong, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
  "Strong reactions", nuclib, opts, True, True)
symmetricFission = REACLIBReactionLibrary(SkyNetRoot
  + "/data/netsu_panov_symmetric_0neut", ReactionType.Strong, False,
  LeptonMode.TreatAllAsDecayExceptLabelEC,
  "Symmetric neutron induced fission with 0 neutrons emitted", nuclib, opts,
  True, True)
spontaneousFission = REACLIBReactionLibrary(SkyNetRoot +
  "/data/netsu_sfis_Roberts2010rates", ReactionType.Strong, False,
  LeptonMode.TreatAllAsDecayExceptLabelEC, "Spontaneous fission", nuclib, opts,
  True, True)

reactionLibraries = [strongReactionLibrary, symmetricFission,
                     spontaneousFission]

screen = SkyNetScreening(nuclib)

net = ReactionNetwork(nuclib, reactionLibraries, helm, screen, opts)
nuclib = net.GetNuclideLibrary()

T = 7
rho = 1.0e9
Ye = 0.4

# calc NSE
nse = NSE(nuclib, helm, screen)
n = nse.CalcFromTemperatureAndDensity(T, rho, Ye)
therm = helm.FromTAndRho(T, rho, n.Y(), nuclib)

dat = np.column_stack((nuclib.Zs(), nuclib.Ns(), n.Y()))
np.savetxt("nse_with_screening_no_DB", dat, fmt="%5i  %5i  %.10e")

print("NSE entropy = %.6e" % therm.S())

Y0 = np.zeros(nuclib.NumNuclides())
Y0[nuclib.NuclideIdsVsNames()["n"]] = 1.0 - Ye
Y0[nuclib.NuclideIdsVsNames()["p"]] = Ye

temperature_vs_time = ConstantFunction(T)
density_vs_time = ConstantFunction(rho)

output = net.Evolve(Y0, 0.0, 1.0e10, temperature_vs_time, density_vs_time,
                    "SkyNet_no_detailed_balance")
