#!/usr/bin/env python

from SkyNet import *
import numpy as np

nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

opts = NetworkOptions()
opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
opts.MassDeviationThreshold = 1.0E-10
opts.IsSelfHeating = True
opts.EnableScreening = True

screen = SkyNetScreening(nuclib)
helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

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

# or use the following code to use FFN rates and weak REACLIB rates pre-computed
# with the <SkyNetRoot>/examples/precompute_reaction_libs.py script

# ffnMesaReactionLibrary = FFNReactionLibrary.ReadFromDisk(
#     "ffnMesa_with_neutrino", opts)
# ffnReactionLibrary = FFNReactionLibrary.ReadFromDisk(
#     "ffn_with_neutrino_ffnMesa", opts)
# weakReactionLibrary = REACLIBReactionLibrary.ReadFromDisk(
#     "weak_REACLIB_with_neutrino_ffnMesa_ffn", opts)

# add neutrino reactions, if desired (will need a neutrino distribution
# function set with net.LoadNeutrinoHistory(...))
# neutrinoLibrary = NeutrinoReactionLibrary(SkyNetRoot
#     + "/data/neutrino_reactions.dat", "Neutrino interactions", nuclib, opts,
#     1.e-2, False, True)

reactionLibraries = [strongReactionLibrary, symmetricFission,
    spontaneousFission, weakReactionLibrary]

# reactionLibraries = [strongReactionLibrary, symmetricFission,
#     spontaneousFission, weakReactionLibrary, neutrinoLibrary,
#     ffnMesaReactionLibrary, ffnReactionLibrary]

net = ReactionNetwork(nuclib, [weakReactionLibrary, strongReactionLibrary,
    symmetricFission, spontaneousFission], helm, screen, opts)

T0 = 6.0    # initial temperature in GK
Ye = 0.01   # initial Ye
s = 10.0    # initial entropy in k_B / baryon
tau = 7.1   # expansion timescale in ms

# run NSE with the temperature and entropy to find the initial density
nse = NSE(net.GetNuclideLibrary(), helm, screen)
nseResult = nse.CalcFromTemperatureAndEntropy(T0, s, Ye)

densityProfile = ExpTMinus3(nseResult.Rho(), tau / 1000.0);

output = net.EvolveSelfHeatingWithInitialTemperature(nseResult.Y(), 0.0, 1.0E9,
    T0, densityProfile, "SkyNet_r-process_python")

NetworkOutput.MakeDatFile("SkyNet_r-process_python.h5")

YvsA = np.array(output.FinalYVsA())
A = np.arange(len(YvsA))

np.savetxt("final_y_r-process_python", np.array([A, YvsA]).transpose(),
    "%6i  %30.20E")
