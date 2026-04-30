#!/usr/bin/env python3
# Parameterized νp-process wind trajectory for vp_process.
# T ~ t^(-2/3), rho ~ t^(-3), r = r0 + v_wind*(t - t_launch),
# qualitatively matching Friedland 2026 fig. 7 (18 Msun GR wind).
# Interim placeholder; replace with real trajectories once available.

import math
import os

Ye = 0.60
t_launch = 1.0 # s
t_end = 30.0   # s
N = 200        # log-spaced points

T0 = 10.0    # GK at t_launch
rho0 = 1.0e8 # g/cm3 at t_launch
r0 = 2.0e6   # cm (PNS surface)
v_wind = 2.0e8 # cm/s (~2000 km/s constant outflow)

times = [t_launch * (t_end / t_launch) ** (k / (N - 1)) for k in range(N)]

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "trajectories")
out_path = os.path.join(out_dir, "vp_wind_param.dat")
os.makedirs(out_dir, exist_ok=True)

with open(out_path, "w") as f:
    f.write("# νp-process parameterized wind trajectory\n")
    f.write("# Ye = {:f}\n".format(Ye))
    f.write("# t_launch = {:.3f} s\n".format(t_launch))
    f.write("# T0 = {} GK, rho0 = {} g/cm3, r0 = {} cm\n".format(T0, rho0, r0))
    f.write("# Columns: time[s]  T[GK]  rho[g/cm3]  r[cm]\n")
    for t in times:
        ratio = t / t_launch
        T = T0 * ratio ** (-2.0 / 3.0)
        rho = rho0 * ratio ** (-3.0)
        r = r0 + v_wind * (t - t_launch)
        f.write("{:.6e}  {:.6e}  {:.6e}  {:.6e}\n".format(t, T, rho, r))

print("Wrote {} rows to {}".format(N, out_path))
