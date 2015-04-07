/// \file NucleiChart.hpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_NUCLEICHART_HPP_
#define SKYNET_MOVIEMAKER_NUCLEICHART_HPP_

#ifdef SKYNET_USE_MOVIE

#include <array>
#include <cstdio>
#include <string>
#include <vector>

#include "MovieMaker/CairoImage.hpp"
#include "MovieMaker/ColorScale.hpp"
#include "MovieMaker/Scale.hpp"
#include "MovieMaker/NucleiChartOptions.hpp"
#include "MovieMaker/MovieData.hpp"

struct Point {
  Point(const double x, const double y):
    x(x),
    y(y) { }

  double x, y;
};

struct PlotCoord {
  PlotCoord(const std::string& str) {
    sscanf(str.c_str(), "%le %le %le %le", &t, &temp, &density, &edot);
  }

  double t;
  double temp;
  double density;
  double edot;
};

class NucleiChart {
public:
  NucleiChart(const NetworkOutput& networkOutput, const double timeStep,
      const ColorScale * const pColorScale, const Scale& abundanceScale,
      const std::string& pathToPlotFiles);

  NucleiChart(const NetworkOutput& networkOutput,
      const double initialTime, const double finalTime, const int numTimes,
      const ColorScale * const pColorScale, const Scale& abundanceScale,
      const std::string& pathToPlotFiles);

  NucleiChart(const MovieData& dat, const ColorScale * const pColorScale,
      const Scale& abundanceScale, const std::string& pathToPlotFiles);

  void DrawFrame(const int timeIdx, const std::string& dir,
      const int saveIdx) const;

  void DrawFrames(const std::string& dir) const;

  static void DrawStandaloneFrame(const std::vector<double>& Y,
      const ThermodynamicState thermoState,
      const NuclideLibrary& nuclideLibrary, const std::string path);

private:
  void MakeBackground(const std::string& pathToPlotFile);

  CairoImage MakeImage() const;

  void DrawHistogram(CairoImage& image, const bool drawBackground,
      const std::vector<double>& binCounts, const int timeIdx) const;

  void DrawLegend(CairoImage& image) const;

  void DrawNumbers(CairoImage& image) const;

  void DrawMissing(CairoImage& image) const;

  // return coordinates of lower left corner of the cell (n,z)
  Point GetPoint(CairoImage& image, const int n, const int z,
      const bool borderCentered) const;

  // return coordinates of center point of the cell (n,z)
  Point GetCellCenterPoint(CairoImage& image, const int n, const int z,
      const bool borderCentered) const;

  void DrawMagicNumberLine(CairoImage& image, const int constNum, const int low,
      const int high, const bool constN) const;

  void DrawMagicNumberAndPeakLines(CairoImage& image) const;

  void DrawConstALine(CairoImage& image, const int a) const;

  void DrawPathLines(CairoImage& image, const int timeIdx) const;

  void DrawStableNuclide(CairoImage& image, const int N, const int Z) const;

  void DrawStableNuclides(CairoImage& image) const;

  void DrawCell(CairoImage& image, const int N, const int Z,
      const std::array<double, 4> fillRGBA,
      const std::array<double, 3> borderRGB, const bool drawBorderOnly) const;

  void DrawNuclide(CairoImage& image, const int N, const int Z,
      const std::array<double, 4> fillRGBA) const {
    DrawCell(image, N, Z, fillRGBA, Opt::borderRGB, false);
  }

  std::vector<double> DrawNuclides(CairoImage& image, const int timeIdx,
      const bool drawBackground) const;

  void DrawPlotDots(CairoImage& image, const int timeIdx) const;

  void DrawText(CairoImage& image, const int timeIdx) const;

  const ColorScale * const mpColorScale;

  const Scale& mrAbundanceScale;

  CairoImage mBackground;
  int mChartWidth;
  int mChartHeight;
  int mWidth;
  int mHeight;
  double mXOffset;
  double mYOffset;
  double mPlotHeight;

  const MovieData mDat;

  std::vector<PlotCoord> mPlotCoords;
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_NUCLEICHART_HPP_
