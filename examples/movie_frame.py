#!/usr/bin/env python

from SkyNet import *
import numpy as np

nuclib = NuclideLibrary.CreateFromWebnucleoXML(SkyNetRoot
    + "/data/webnucleo_nuc_v2.0.xml")

screen = SkyNetScreening(nuclib)
helm = HelmholtzEOS(SkyNetRoot + "/data/helm_table.dat")

T0 = 7 # GK
s = 10 # k_B / baryon
Ye = 0.1

nse = NSE(nuclib, helm, screen)
nseResult = nse.CalcFromTemperatureAndEntropy(T0, s, Ye)
y = np.array(nseResult.Y())

therm = helm.FromTAndRho(T0, nseResult.Rho(), y, nuclib)

NucleiChart.DrawStandaloneFrame(y, therm, nuclib, "NSE_7GK_10kB_Ye_0.1.png")
