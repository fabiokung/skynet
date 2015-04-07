/// \file movie.cpp
/// \author jlippuner
/// \since Oct 15, 2014
///
/// \brief
///
///

#include <string>

#include "MovieMaker/ColorScales.hpp"
#include "MovieMaker/NucleiChart.hpp"
#include "Network/NetworkOutput.hpp"

int main(int argc, char** argv) {
  if (argc != 2) {
    printf("movie expects one argument (particleXXX)\n");
    return EXIT_FAILURE;
  }

  std::string h5Path(argv[1]);

  // TODO make this more robust
  auto lastSlash = h5Path.find_last_of("/");

  std::string rootPath = ".";
  if (lastSlash != std::string::npos)
    rootPath = h5Path.substr(0, lastSlash + 1);

  ColorScale * pColorScale = new WhiteBlueGreenYellowRedColorScale(
  //ColorScale * pColorScale = new WhiteBlackColorScale(
      1.0E-12, 1.0E-2,
      true, { 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2 },
      { "10^{-12}", "10^{-11}", "10^{-10}", "10^{-9}", "10^{-8}", "10^{-7}", "10^{-6}", "10^{-5}",
          "10^{-4}", "10^{-3}", "10^{-2}" });

  Scale abundanceScale(1E-12, 1E-2, true, { 1E-12, 1E-10, 1E-8, 1E-6, 1E-4, 1E-2 },
      { "10^{-12}", "10^{-10}", "10^{-8}", "10^{-6}", "10^{-4}", "10^{-2}" });

  auto out = NetworkOutput::ReadFromFile(h5Path);

  double timeStep = 0.02;
  MovieData movieData(out, std::max(1.0E-3, out.Times()[0]), out.Times()[out.Times().size() - 1],
      timeStep);
  NucleiChart chart(movieData, pColorScale, abundanceScale, rootPath);

//  chart.DrawFrame(0, rootPath, 0);
//  chart.DrawFrame(1, rootPath, 1);
//  chart.DrawFrame(100, rootPath, 100);
//  chart.DrawFrame(200, rootPath, 200);
//  chart.DrawFrame(584, rootPath, 584);
//  chart.DrawFrame(1502, rootPath, 1502);

  chart.DrawFrames(rootPath + "/chart_frames/");

  delete pColorScale;

  return EXIT_SUCCESS;
}

