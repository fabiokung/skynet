#!/usr/bin/env python

from SkyNet import *
import numpy as np
import time

nucs = open("nuclide_names", "r").readlines()
nucs = [n.strip() for n in nucs]

nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml", nucs)

ni56 = nuclib.NuclideIdsVsNames()["ni56"]
co55 = nuclib.NuclideIdsVsNames()["co55"]
fe55 = nuclib.NuclideIdsVsNames()["fe55"]
p = nuclib.NuclideIdsVsNames()["p"]

helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

scr = SkyNetScreening(nuclib)

Ye = 0.5   # initial Ye
rho = 5.0e8

nse_opts = NSEOptions()
nse_opts.SetDoScreening(False)

nse = NSE(nuclib, helm, None, nse_opts)
nse_scr = NSE(nuclib, helm, scr)

Ts = np.linspace(3, 12, 300)

print("# [1] = Temperature [GK]")
print("# [2] = Ni-56 mass fraction, no screening")
print("# [3] = Co-55 mass fraction, no screening")
print("# [4] = Fe-55 mass fraction, no screening")
print("# [5] = proton mass fraction, no screening")
print("# [6] = Ni-56 mass fraction, screening")
print("# [7] = Co-55 mass fraction, screening")
print("# [8] = Fe-55 mass fraction, screening")
print("# [9] = proton mass fraction, screening")

for t in Ts:
  nsc = nse.CalcFromTemperatureAndDensity(t, rho, Ye).Y()
  sc = nse_scr.CalcFromTemperatureAndDensity(t, rho, Ye).Y()

  screen = scr.Compute(helm.FromTAndRho(t, rho, sc, nuclib), sc, nuclib)

  print("%5.2f %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e" % (
    t, nsc[ni56] * 56., nsc[co55] * 55., nsc[fe55] * 55., nsc[p],
       sc[ni56] * 56., sc[co55] * 55., sc[fe55] * 55., sc[p], screen.Mus()[p]))
