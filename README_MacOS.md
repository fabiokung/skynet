# Building SkyNet on macOS (Apple Silicon)

Everything builds from Homebrew and the Xcode Command Line Tools. The full
example networks (~8000 nuclides) need a sparse solver; on Apple Silicon that is
Pardiso (free for academic use). The `macos` preset uses dense Armadillo
instead — no license, but only practical for small networks.

Targets Apple Silicon (Apple Clang, Homebrew Python 3 / SWIG 4).

## Packages

```bash
xcode-select --install   # if the Command Line Tools aren't installed yet
brew install cmake gsl hdf5 boost swig python libomp armadillo \
  cairomm@1.14 libsigc++@2 ffmpeg
```

`cairomm@1.14` and `libsigc++@2` are the keg-only versions the movie code builds
against (the cairomm-1.0 / sigc++-2.0 APIs).

## Pardiso

Pardiso is a separate download with a free academic license:

1. Register at panua.org and download the Apple Silicon `libpardiso.dylib`.
2. Put the library at `~/lib/pardiso/libpardiso.dylib` and the license at
   `~/panua.lic`.

The dylib has a bare install name, so `DYLD_LIBRARY_PATH` must point at it when
running binaries directly (the `pardiso` test preset sets this for `ctest`).

## Build

```bash
PKG_CONFIG_PATH="/opt/homebrew/opt/cairomm@1.14/lib/pkgconfig:/opt/homebrew/opt/libsigc++@2/lib/pkgconfig" \
  cmake --preset pardiso -DENABLE_MOVIE=ON -DREQUIRE_MOVIE=ON
cmake --build build -j$(sysctl -n hw.ncpu)
cmake --build build --target install
```

`--preset pardiso` enables the SWIG bindings; the `PKG_CONFIG_PATH` +
`ENABLE_MOVIE` add the Cairo movie tools. Install is required before running:
`SkyNetRoot` points at the install tree, so binaries read nuclide/EOS data from
`skynet_install/data`.

For a license-free build without sparse solving, use `cmake --preset macos`
(Armadillo + bindings, small networks only).

## Run a network

`r-process` is self-contained and writes `SkyNet_r-process.h5`. Pardiso needs
`DYLD_LIBRARY_PATH` at runtime:

```bash
mkdir run && cd run
DYLD_LIBRARY_PATH="$HOME/lib/pardiso" ../skynet_install/examples/r-process
```

It loads the full network and self-heats to t = 1e9 s, so it runs a while; the
h5 is readable mid-run.

## Make a movie

The chart movie needs no Python. The optional T/ρ/heating overlay needs
matplotlib, scipy, and pycairo — Homebrew's Python is externally managed, so put
them in a venv:

```bash
python3 -m venv .venv && .venv/bin/pip install matplotlib scipy pycairo
```

From the directory holding the h5:

```bash
export DYLD_LIBRARY_PATH="$HOME/lib/pardiso"
export PYTHONPATH=$PWD/../skynet_install/lib
.venv/bin/python ../skynet_install/examples/plot_temp_rho_edot.py SkyNet_r-process.h5
mkdir chart_frames
../skynet_install/examples/movie SkyNet_r-process.h5
ffmpeg -y -c:v png -i "chart_frames/chart_%06d.png" \
  -pix_fmt yuv420p -c:v libx264 -preset slow -crf 1 SkyNet_r-process.mp4
```

Example scripts carry a `#!/usr/bin/env python` shebang; run them through an
explicit interpreter (`.venv/bin/python <script>`). The overlay is optional —
`movie` renders the chart without `plot_background.png`.
```
