# Building SkyNet on Ubuntu LTS

Everything builds from apt. The `ubuntu` preset uses the Intel oneAPI MKL sparse
solver and enables the SWIG bindings and movie tools. MKL is free and fast on the
full ~8000-nuclide example networks, but x86_64-only — on ARM, or to skip the
Intel apt repo, use the Trilinos KLU build at the bottom (slower).

Targets Ubuntu 24.04 (cmake 3.28, SWIG 4.2, Python 3.12, gfortran 13).

## Packages

```bash
sudo apt-get install -y \
  build-essential cmake gfortran git wget gpg ca-certificates \
  libhdf5-dev libgsl-dev libboost-filesystem-dev libboost-serialization-dev \
  swig python3-dev python3-numpy \
  libfreetype-dev libcairo2-dev libcairomm-1.0-dev libsigc++-2.0-dev ffmpeg \
  python3-matplotlib python3-scipy python3-cairo   # last line: movies only
```

Intel MKL (installs to `/opt/intel/oneapi/mkl/latest`, where the preset looks):

```bash
wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | sudo gpg --dearmor -o /usr/share/keyrings/oneapi-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  | sudo tee /etc/apt/sources.list.d/oneAPI.list
sudo apt-get update && sudo apt-get install -y intel-oneapi-mkl-devel
```

MKL is linked statically, so there's no `LD_LIBRARY_PATH` at runtime.

## Build

```bash
cmake --preset ubuntu
cmake --build build -j$(nproc)
cmake --build build --target install
```

Install is required before running: `SkyNetRoot` points at the install tree, so
binaries read nuclide/EOS data from `skynet_install/data`. Re-install only when
`data/` changes.

## Run + movie

`r-process` is self-contained (no input files) and writes `SkyNet_r-process.h5`.
It loads the full network and self-heats to t = 1e9 s, so it runs a while; the h5
is readable mid-run.

```bash
mkdir run && cd run
../skynet_install/examples/r-process

export PYTHONPATH=$PWD/../skynet_install/lib
python3 ../skynet_install/examples/plot_temp_rho_edot.py SkyNet_r-process.h5
mkdir chart_frames
../skynet_install/examples/movie SkyNet_r-process.h5
ffmpeg -y -c:v png -i "chart_frames/chart_%06d.png" \
  -pix_fmt yuv420p -c:v libx264 -preset slow -crf 1 SkyNet_r-process.mp4
```

Example scripts carry a `#!/usr/bin/env python` shebang; run them as `python3
<script>` (Ubuntu has no `python`). The overlay is optional — `movie` renders
without `plot_background.png`.

## Trilinos KLU (license-free, any arch)

Replaces MKL. SkyNet's `FindTrilinos` wants every advertised component, so install
the whole set:

```bash
sudo apt-get install -y 'libtrilinos-*-dev' libsuitesparse-dev
cmake --preset ubuntu-trilinos
cmake --build build -j$(nproc)
cmake --build build --target install
```
