#!/usr/bin/env python

from SkyNet import *

nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

reaclib = REACLIB(SkyNetRoot + "/data/reaclib")

opts = NetworkOptions()
opts.ConvergenceCriterion = NetworkConvergenceCriterion.Mass
opts.MassDeviationThreshold = 1.0E-10
opts.IsSelfHeating = True
opts.CalculateStrongInverseRates = True
opts.DisableStdoutOutput = True

out = NetworkOutput.CreateNew("precompute_reaction_libs", nuclib, False)

ffnMesaReactionLibrary = FFNReactionLibrary(
    SkyNetRoot + "/data/FFN_mesa_weak_rates.h5",
    ReactionType.Weak, "Weak MESA reactions", nuclib, opts, True)
print("Loaded FFN MESA library")

ffnMesaReactionLibrary.Dump("ffnMesa")
print("Dumped FFN MESA library")

ffnReactionLibrary = FFNReactionLibrary(
    SkyNetRoot + "/data/FFN_full_weak_rates.h5",
    ReactionType.Weak, "Weak FFN reactions", nuclib, opts, True)
print("Loaded FFN library")

ffnReactionLibrary.RemoveReactions(ffnMesaReactionLibrary.Reactions())
print("Trimmed FFN library")

ffnReactionLibrary.Dump("ffn_with_ffnMesa")
print("Dumped FFN library")

weakReactionLibrary = REACLIBReactionLibrary(SkyNetRoot + "/data/reaclib",
    ReactionType.Weak, False, LeptonMode.TreatAllAsDecayExceptLabelEC,
    "Weak REACLIB reactions", nuclib, opts, True)
print("Loaded weak library")

weakReactionLibrary.RemoveReactions(ffnMesaReactionLibrary.Reactions())
weakReactionLibrary.RemoveReactions(ffnReactionLibrary.Reactions())
print("Trimmed weak library")

weakReactionLibrary.Dump("weak_REACLIB_with_ffnMesa_ffn")
print("Dumped weak library")

ffnMesaReactionLibrary.PrintRates(out)
ffnReactionLibrary.PrintRates(out)
weakReactionLibrary.PrintRates(out)
