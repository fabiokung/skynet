# MKL build (Linux)

The `mkl` CMake preset links the Intel oneAPI MKL sparse solver. Linux-only
(hidden on macOS by a preset condition); on Apple Silicon use the `pardiso`
preset. Both give the same Y(A) to round-off (~6e-7). For a from-scratch Ubuntu
setup (apt packages, bindings, movies), see README_Ubuntu.md.

## Prerequisites

- Intel oneAPI MKL at `/opt/intel/oneapi/mkl/latest` (set `MKLROOT` if elsewhere)
- gcc/gfortran — no Intel compiler needed
- MKL is linked statically, so there's no `LD_LIBRARY_PATH` at runtime

## Build

```bash
cmake --preset mkl
cmake --build build -j$(nproc)
cmake --build build --target install   # populates skynet_install/data
```

Install before running: `SkyNetRoot` points at the install tree, so binaries read
nuclide/EOS data from `skynet_install/data`.

## Multithreaded build

`FindMKL.cmake` exposes `MKL_THREADING` (default `sequential`). Set it to `gnu` to
link `mkl_gnu_thread` + libgomp and thread the sparse factorization. Use a separate
build dir so the sequential build stays intact:

```bash
MKLROOT=/opt/intel/oneapi/mkl/latest cmake -S . -B build_mt \
  -DSKYNET_MATRIX_SOLVER=mkl -DMKL_THREADING=gnu \
  -DCMAKE_INSTALL_PREFIX=$PWD/skynet_install
cmake --build build_mt -j$(nproc)
```

Thread count is set at runtime via `OMP_NUM_THREADS`. Only the factorization is
threaded — the Jacobian assembly and rate evaluation are single-threaded — so the
full-run speedup is solver-bound and plateaus by ~8 threads (roughly 1.5–1.7x).
Running two cases at 8 threads each beats one 16-thread run on aggregate throughput.
