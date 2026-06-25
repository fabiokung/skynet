#!/bin/bash
# Build SkyNet's Cairo MovieMaker and render a nuclide-chart movie for one run.
# Follows README_make_movie: overlay plot -> chart frames -> ffmpeg.
#
# Needs brew cairomm@1.14 (= the cairomm-1.0 API series, so the code builds
# unchanged), libsigc++@2, ffmpeg. The movie binary lives in a separate
# build-movie/ tree (ENABLE_MOVIE=ON); the main build/ has movies off.
#
# Each run renders in its own <stem>_movie/ dir because NucleiChart reads
# plot_coords.txt / plot_background.png and writes chart_frames/ next to the
# h5 -- two runs in one dir would clobber each other.
#
# Usage: ./make_movie.sh <run.h5> [output.mp4]
#   OVERLAY=0 ./make_movie.sh ...   skip the T/rho/heating-rate overlay
set -e

cd "$(dirname "$0")"

CAIROMM=$(brew --prefix cairomm@1.14)
SIGCPP=$(brew --prefix libsigc++@2)
BREW=$(brew --prefix)
export PKG_CONFIG_PATH="$CAIROMM/lib/pkgconfig:$SIGCPP/lib/pkgconfig:$BREW/lib/pkgconfig:$BREW/share/pkgconfig"
export DYLD_LIBRARY_PATH="$HOME/lib/pardiso${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"

# movie-enabled build (configure once, then incremental)
[ -f ../../build-movie/CMakeCache.txt ] || \
  cmake --preset pardiso -B ../../build-movie -DENABLE_MOVIE=ON -DREQUIRE_MOVIE=ON
cmake --build ../../build-movie --target movie

h5=$1
[ -n "$h5" ] || { echo "usage: $0 <run.h5> [output.mp4]" >&2; exit 1; }
mp4=${2:-${h5%.h5}.mp4}
name=$(basename "$h5")
work=$(dirname "$h5")/$(basename "${h5%.h5}")_movie
rm -rf "$work"
mkdir -p "$work/chart_frames"
ln -sf "../$name" "$work/$name"

# T/rho/heating-rate background + time-marker coords (NetworkOutput only, so the
# non-movie build/ bindings are enough). Needs scipy + pycairo in the venv.
if [ "${OVERLAY:-1}" != "0" ]; then
  cache=${TMPDIR:-/tmp}/skynet_movie_cache
  mkdir -p "$cache/fontconfig"
  MPLCONFIGDIR="$cache" XDG_CACHE_HOME="$cache" \
    PYTHONPATH="$PWD/../../build" \
    .venv/bin/python ../plot_temp_rho_edot.py "$work/$name"
fi

../../build-movie/examples/movie/movie "$work/$name"

ffmpeg -y -c:v png -i "$work/chart_frames/chart_%06d.png" \
  -pix_fmt yuv420p -c:v libx264 -preset slow -crf 1 "$mp4"
echo "wrote $mp4"
