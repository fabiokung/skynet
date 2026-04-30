/// \file NucleiChart.cpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#include "MovieMaker/NucleiChart.hpp"

#include <cassert>
#include <cmath>
#include <fstream>

#include "BuildInfo.hpp"
#include "MovieMaker/ColorScales.hpp"

NucleiChart::NucleiChart(const NetworkOutput& networkOutput,
    const double timeStep, const ColorScale * const pColorScale,
    const Scale& abundanceScale, const std::string& pathToPlotFiles) :
    NucleiChart(MovieData(networkOutput, timeStep), pColorScale, abundanceScale,
        pathToPlotFiles) {}

NucleiChart::NucleiChart(const NetworkOutput& networkOutput,
    const double initialTime, const double finalTime, const int numTimes,
    const ColorScale * const pColorScale, const Scale& abundanceScale,
    const std::string& pathToPlotFiles) :
    NucleiChart(MovieData(networkOutput, initialTime, finalTime, numTimes),
        pColorScale, abundanceScale, pathToPlotFiles) {
}

NucleiChart::NucleiChart(const MovieData& dat,
    const ColorScale * const pColorScale, const Scale& abundanceScale,
    const std::string& pathToPlotFile):
    mpColorScale(pColorScale),
    mrAbundanceScale(abundanceScale),
    mBackground(0, 0),
    mDat(dat) {

  assert(Opt::cellSize > 0);
  assert(Opt::lineWidth >= 0);

  // read plot coords
  if (pathToPlotFile != "") {
    std::ifstream istm(pathToPlotFile + "/plot_coords.txt");
    std::string line;
    while (std::getline(istm, line)) {
      if (line[0] == '#')
        continue;
      mPlotCoords.push_back(PlotCoord(line));
    }
  }

  // make background
  MakeBackground(pathToPlotFile);

  //mBackground.Save("background.png");
}

void NucleiChart::DrawFrame(const int timeIdx, const std::string& dir,
    const int saveIdx) const {
  if (timeIdx < mDat.FirstTimeIndex()) {
    return;
  }

  auto image = MakeImage();
  char buf[128];

  std::vector<double> binCounts = DrawNuclides(image, timeIdx, false);
  DrawMagicNumberAndPeakLines(image);
  DrawPathLines(image, timeIdx);
  DrawHistogram(image, false, binCounts, timeIdx);
  DrawPlotDots(image, timeIdx);
  DrawText(image, timeIdx);

  sprintf(buf, "chart_%06d.png", saveIdx);
  image.Save(dir + "/" + std::string(buf));
  printf("done frame %06d\n", saveIdx);
}

void NucleiChart::DrawFrames(const std::string& dir) const {
  for (int t = 0; t < mDat.NumTimes(); ++t) {
    DrawFrame(t, dir, t - mDat.FirstTimeIndex());
  }

  int numFrames = mDat.NumTimes() - mDat.FirstTimeIndex();
  // round up to multiple of 25
  int numFrames25 = 25 * (numFrames / 25 + 1);

  // pad with last frame until we have a multiple of 25
  for (int i = 0; i < (numFrames25 - numFrames); ++i) {
    DrawFrame(mDat.NumTimes() - 1, dir,
        mDat.NumTimes() + i - mDat.FirstTimeIndex());
  }
}

void NucleiChart::DrawStandaloneFrame(const std::vector<double>& Y,
    const ThermodynamicState thermoState,
    const NuclideLibrary& nuclideLibrary, const std::string path) {
  auto out = NetworkOutput::CreateNew("/tmp/skynet_output", nuclideLibrary,
      true);

  out.AddEntry( 0.1, { 1.0, "init" }, 0, 0, thermoState, 0.0, Y, 1.0);
  out.AddEntry(10.0, { 1.0, "init" }, 0, 0, thermoState, 0.0, Y, 1.0);

  ColorScale * pColorScale = new WhiteBlueGreenYellowRedColorScale(
  //ColorScale * pColorScale = new WhiteBlackColorScale(
      1.0E-12, 1.0E-2,
      true, { 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5, 1E-4, 1E-3, 1E-2 },
      { "10^{-12}", "10^{-11}", "10^{-10}", "10^{-9}", "10^{-8}", "10^{-7}", "10^{-6}", "10^{-5}",
          "10^{-4}", "10^{-3}", "10^{-2}" });

  Scale abundanceScale(1E-12, 1E-2, true, { 1E-12, 1E-10, 1E-8, 1E-6, 1E-4, 1E-2 },
      { "10^{-12}", "10^{-10}", "10^{-8}", "10^{-6}", "10^{-4}", "10^{-2}" });

  NucleiChart chart(out, 1.0, 2.0, 2, pColorScale, abundanceScale, "");

  auto image = chart.MakeImage();
  int timeIdx = 0;

  std::vector<double> binCounts = chart.DrawNuclides(image, timeIdx, false);
  chart.DrawMagicNumberAndPeakLines(image);
  chart.DrawPathLines(image, timeIdx);
  chart.DrawHistogram(image, false, binCounts, timeIdx);
  chart.DrawPlotDots(image, timeIdx);
  chart.DrawText(image, timeIdx);

  image.Save(path);
}

void NucleiChart::MakeBackground(const std::string& pathToPlotFile) {
  double totalBorder = Opt::lowerLeftBorder + Opt::upperRightBorder;

  // neutron axis goes from N = 0 to N = mMaxN
  mChartWidth = (Opt::cellSize + Opt::lineWidth / 2.0) * (mDat.MaxN() + 1)
      + Opt::lineWidth / 2.0;
  int width = mChartWidth
      + int((Opt::cellSize + Opt::lineWidth / 2.0) * totalBorder);

  // proton axis goes from Z = 0 to Z = mMaxZ
  mChartHeight = (Opt::cellSize + Opt::lineWidth / 2.0) * (mDat.MaxZ() + 1)
      + Opt::lineWidth / 2.0;
  int height = mChartHeight
      + int((Opt::cellSize + Opt::lineWidth / 2.0) * totalBorder);

  // to encode the video with some advanced codecs (e.g. h264), the width and
  // heights must be multiples of 16
  mWidth = 16 * (width / 16 + 1);
  mHeight = 16 * (height / 16 + 1);
  int dx = double(mWidth - width) / 2.0;
  int dy = double(mHeight - height) / 2.0;
  mXOffset = dx + (Opt::cellSize + Opt::lineWidth / 2.0) * Opt::lowerLeftBorder;
  mYOffset = -dy
      - (Opt::cellSize + Opt::lineWidth / 2.0) * Opt::lowerLeftBorder;

  mBackground = CairoImage(mWidth, mHeight, Opt::backgroundRGB);
  mBackground.Context()->set_line_width(Opt::lineWidth);

  // load plot background (optional — skip silently if the file doesn't exist)
  mPlotHeight = 0.0;
  if (pathToPlotFile != "") {
    std::ifstream plotBgCheck(pathToPlotFile + "/plot_background.png");
    if (plotBgCheck.good()) {
      plotBgCheck.close();
      CairoImage plotBackground(0, 0);
      plotBackground.Surface() = Cairo::ImageSurface::create_from_png(
          pathToPlotFile + "/plot_background.png");
      mPlotHeight = plotBackground.Height();
      mBackground.Context()->save();
      mBackground.Context()->translate(mXOffset / 2.0,
          mHeight - mChartHeight + mYOffset - 5.0);
      mBackground.Context()->set_source(plotBackground.Surface(), 0.0, 0.0);
      mBackground.Context()->paint();
      mBackground.Context()->restore();
    }
  }

  // load logo
  CairoImage logo(0, 0);
  logo.Surface() = Cairo::ImageSurface::create_from_png(
      SkyNetRoot + "/data/logo/SkyNet_logo.png");
  mBackground.Context()->save();
  mBackground.Context()->translate(mXOffset / 2.0,
        mHeight - mChartHeight + mYOffset -5.0 + mPlotHeight);
  mBackground.Context()->set_source(logo.Surface(), 0.0, 0.0);
  mBackground.Context()->paint();
  mBackground.Context()->restore();

  mBackground.Context()->translate(mXOffset, mYOffset);

  DrawNuclides(mBackground, 0, true);
  DrawLegend(mBackground);
  DrawNumbers(mBackground);
  DrawHistogram(mBackground, true, {}, 0);

  // add logo
  // TODO improve this, add actual logo
  mBackground.Context()->set_font_size(0.85 * Opt::fontSize);
  double logoX = 250.0;
  mBackground.DrawText(logoX, mHeight - 1.0, "Made with SkyNet by Jonas Lippuner",
        CairoImage::HorizontalTextPosition::Left);
}

CairoImage NucleiChart::MakeImage() const {
  CairoImage image(mWidth, mHeight, Opt::backgroundRGB);
  image.Context()->set_source(mBackground.Surface(), 0.0, 0.0);
  image.Context()->paint();
  image.Context()->set_line_width(Opt::lineWidth);

  image.Context()->translate(mXOffset, mYOffset);

  return image;
}

void NucleiChart::DrawHistogram(CairoImage& image, const bool drawBackground,
    const std::vector<double>& binCounts, const int timeIdx) const {
  if (!drawBackground) {
    assert(binCounts.size() == Opt::leftBinEdges.size());
    assert(Opt::binNames.size() == Opt::leftBinEdges.size());
  }

  // TODO
  double x0 = 570.0;
  double y0 = image.Height() - 1.0 - 30.0;
  double totalHeight = 150.0;

  image.Context()->save();
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  image.Context()->set_line_width(2.0);

  Cairo::RefPtr<Cairo::ToyFontFace> font;
  if (drawBackground) {
    font = Cairo::ToyFontFace::create(Opt::fontFamily, Cairo::FONT_SLANT_NORMAL,
        Cairo::FONT_WEIGHT_NORMAL);
    image.Context()->set_font_face(font);
    image.Context()->set_font_size(Opt::fontSize);
  } else {
    image.Context()->set_antialias(Cairo::ANTIALIAS_NONE);
  }

  // draw bins
  for (unsigned int i = 0; i < Opt::leftBinEdges.size(); ++i) {
    double x = x0 + double(Opt::leftBinEdges[i] - 1) * Opt::binWdith;

    int rightEdge = (i == Opt::leftBinEdges.size() - 1) ? mDat.MaxA() :
        Opt::leftBinEdges[i+1];
    int widthInA = rightEdge - Opt::leftBinEdges[i];

    double width = double(widthInA) * Opt::binWdith;

    if (drawBackground) {
      image.DrawText(x + width, y0 + Opt::fontSize / 2.0,
          std::to_string(rightEdge), CairoImage::HorizontalTextPosition::Center,
          CairoImage::VerticalTextPosition::Top);

      image.Context()->set_font_size(Opt::fontSize * 0.85);
      image.DrawText(x + width / 2.0, y0 + Opt::fontSize * 1.8,
          Opt::binNames[i], CairoImage::HorizontalTextPosition::Center,
          CairoImage::VerticalTextPosition::Top);
      image.Context()->set_font_size(Opt::fontSize);
    } else {
      double height = totalHeight * binCounts[i];

      std::array<double, 4> rgba = Opt::binColors[i];
      image.Context()->set_source_rgba(rgba[0], rgba[1], rgba[2], rgba[3]);
      image.Context()->rectangle(x + 1.0, y0 - 1.0, width - 2.0, -height + 2.0);
      image.Context()->fill();

      image.Context()->set_source_rgb(0.0, 0.0, 0.0);
      image.Context()->rectangle(x, y0, width, -height);
      image.Context()->stroke();
    }
  }

  // draw abundance
  if (!drawBackground) {
    image.Context()->save();
    image.Context()->set_source_rgb(0.0, 0.0, 1.0);
    image.Context()->set_line_width(Opt::lineWidth / 2.0);
    image.Context()->set_antialias(Cairo::ANTIALIAS_DEFAULT);

    double abund = mDat.AbundanceVsAVsTLinInterp()[1][timeIdx];
    image.Context()->move_to(x0, y0 - mrAbundanceScale(abund) * totalHeight);

    for (int a = 2; a <= mDat.MaxA(); ++a) {
      double abund = mDat.AbundanceVsAVsTLinInterp()[a][timeIdx];
      image.Context()->line_to(x0 + double(a - 1) * Opt::binWdith,
          y0 - mrAbundanceScale(abund) * totalHeight);
    }

    image.Context()->stroke();
    image.Context()->restore();
  }

  // draw y-axis
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);

  image.Context()->move_to(x0 - 1.0, y0 + 1.0 - totalHeight);
  image.Context()->line_to(x0 - 1.0, y0 + 1.0);

  image.Context()->move_to(x0 - 1.0, y0);
  image.Context()->line_to(x0 + double(mDat.MaxA()) * Opt::binWdith - 1.0,
      y0);

  image.Context()->move_to(x0 + double(mDat.MaxA()) * Opt::binWdith - 1.0,
      y0 + 1.0);
  image.Context()->line_to(x0 + double(mDat.MaxA()) * Opt::binWdith - 1.0,
      y0 + 1.0 - totalHeight);

  image.Context()->stroke();

  if (drawBackground) {
    image.DrawText(x0 - Opt::fontSize / 2.0, y0, "0",
        CairoImage::HorizontalTextPosition::Right,
        CairoImage::VerticalTextPosition::Bottom);
    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - totalHeight, "1",
        CairoImage::HorizontalTextPosition::Right);

//    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - totalHeight / 2.0, "0.5",
//        CairoImage::HorizontalTextPosition::Right);
    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - 0.2 * totalHeight, "0.2",
        CairoImage::HorizontalTextPosition::Right);
    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - 0.4 * totalHeight, "0.4",
        CairoImage::HorizontalTextPosition::Right);
    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - 0.6 * totalHeight, "0.6",
        CairoImage::HorizontalTextPosition::Right);
    image.DrawText(x0 - Opt::fontSize / 2.0, y0 - 0.8 * totalHeight, "0.8",
        CairoImage::HorizontalTextPosition::Right);

    for (unsigned int i = 0; i < mrAbundanceScale.Tics().size(); ++i) {
      double tic = mrAbundanceScale.Tics()[i];
      double yTic = mrAbundanceScale.UseLog() ? log(tic) : tic;
      double y = (yTic - mrAbundanceScale.XMin()) / mrAbundanceScale.Dx()
          * (-totalHeight);

      image.DrawText(x0 + double(mDat.MaxA()) * Opt::binWdith - 1.0 + 3.0,
          y0 + y, mrAbundanceScale.TicLabels()[i],
          CairoImage::HorizontalTextPosition::Left,
          i == 0 ? CairoImage::VerticalTextPosition::Bottom :
              CairoImage::VerticalTextPosition::Middle);
    }

    font = Cairo::ToyFontFace::create(Opt::fontFamily, Cairo::FONT_SLANT_NORMAL,
        Cairo::FONT_WEIGHT_BOLD);
    image.Context()->set_font_face(font);
    image.DrawText(x0 + double(mDat.MaxA() + 10) * Opt::binWdith,
        y0 + Opt::fontSize * 1.5, "A",
        CairoImage::HorizontalTextPosition::Center,
        CairoImage::VerticalTextPosition::Top);

    image.Context()->save();
    image.Context()->translate(x0 - Opt::fontSize * 3.0,
        y0 - totalHeight / 2.0);
    image.Context()->rotate_degrees(-90.0);
    image.DrawText(0.0, 0.0, "Mass fraction");
    image.Context()->restore();

    image.Context()->save();
    image.Context()->translate(
        x0 + double(mDat.MaxA()) * Opt::binWdith + Opt::fontSize * 3.5,
        y0 - totalHeight / 2.0);
    image.Context()->rotate_degrees(-90.0);
    image.Context()->set_source_rgb(0.0, 0.0, 1.0);
    image.DrawText(0.0, 0.0, "Abundance");
    image.Context()->restore();
  }

  image.Context()->restore();
}

void NucleiChart::DrawLegend(CairoImage& image) const {
  //Point pLowerLeft = GetPoint(image, 196, 20, false);
  //Point pUpperRight = GetPoint(image, 199, 70, false);

  double right = mChartWidth - Opt::colorScaleLeftBorder;
  double left = right - Opt::colorScaleWidth;
  //double bottom = mHeight - mChartHeight + mYOffset;
  double bottom = mHeight;

  int nForTop = std::max(mDat.MaxN() - 20, 0);
  double top = GetPoint(image, 0, mDat.ZLimitsVsN()[nForTop].min, false).y;
  if (nForTop > 0)
    top = std::min(top, GetPoint(image, 0, mDat.ZLimitsVsN()[nForTop - 1].min,
        false).y);
  if (nForTop > 1)
    top = std::min(top, GetPoint(image, 0, mDat.ZLimitsVsN()[nForTop - 2].min,
        false).y);

  top += Opt::colorScaleTopBorder;

  Point pLowerLeft(left, bottom);
  Point pUpperRight(right, top);
  double width = pUpperRight.x - pLowerLeft.x;
  double height = pUpperRight.y - pLowerLeft.y;

  // make color scale
  int steps = int(-height);
  for (int i = 0; i < steps; ++i) {
    double x = mpColorScale->GetScale().XMin() + double(i) / (-height) *
        mpColorScale->GetScale().Dx();
    if (mpColorScale->GetScale().UseLog())
      x = exp(x);

    std::array<double, 4> rgba = mpColorScale->GetRGBA(x);
    image.Context()->set_source_rgb(rgba[0], rgba[1], rgba[2]);
    image.Context()->rectangle(pLowerLeft.x, pLowerLeft.y - double(i+1),
        width, 1.0);
    image.Context()->fill();
  }

  // make tics
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);

  Cairo::RefPtr<Cairo::ToyFontFace> font = Cairo::ToyFontFace::create(
      Opt::fontFamily, Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
  image.Context()->set_font_face(font);
  image.Context()->set_font_size(Opt::fontSize);

  for (unsigned int i = 0; i < mpColorScale->GetScale().Tics().size(); ++i) {
    double tic = mpColorScale->GetScale().Tics()[i];
    double yTic = mpColorScale->GetScale().UseLog() ? log(tic) : tic;
    double y = (yTic - mpColorScale->GetScale().XMin())
        / mpColorScale->GetScale().Dx() * (-height);

    image.Context()->move_to(pLowerLeft.x, pLowerLeft.y - y);
    image.Context()->line_to(pLowerLeft.x + width / 4.0, pLowerLeft.y - y);
    image.Context()->move_to(pUpperRight.x, pLowerLeft.y - y);
    image.Context()->line_to(pUpperRight.x - width / 4.0, pLowerLeft.y - y);
    image.Context()->stroke();

    image.DrawText(pLowerLeft.x - width / 2.0, pLowerLeft.y - y,
        mpColorScale->GetScale().TicLabels()[i],
        CairoImage::HorizontalTextPosition::Right,
        CairoImage::VerticalTextPosition::Middle);
  }

  // make border
  image.Context()->rectangle(pLowerLeft.x, pLowerLeft.y, width, height);
  image.Context()->stroke();

  // make legend
  // TODO
  int n = 180;
  int z = 65;
  int dz = 5;

  image.Context()->save();
  image.Context()->set_line_width(Opt::lineWidth / 2.0);
  image.Context()->set_antialias(Cairo::ANTIALIAS_NONE);

  DrawNuclide(image, n, z, mpColorScale->GetRGBARange01(0.0));
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  Point p = GetPoint(image, n + 2, z, false);
  image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5, "Unstable nuclide",
      CairoImage::HorizontalTextPosition::Left);

  z -= dz;
  DrawStableNuclide(image, n, z);
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  p = GetPoint(image, n + 2, z, false);
  image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5, "Stable nuclide",
      CairoImage::HorizontalTextPosition::Left);

  z -= dz;
  DrawNuclide(image, n, z, Opt::missingFillRGBA);
  p = GetPoint(image, n + 2, z, false);
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5, "Missing nuclide",
      CairoImage::HorizontalTextPosition::Left);

  image.Context()->restore();

  z -= dz;
  Point p0 = GetPoint(image, n, z, true);
  p0.y += Opt::lineWidth / 2.0;
  Point p1 = GetPoint(image, n+1, z+1, true);

  image.Context()->move_to(p0.x, p0.y);
  image.Context()->line_to(p0.x, p1.y);
  image.Context()->stroke();

  image.Context()->move_to(p1.x, p0.y);
  image.Context()->line_to(p1.x, p1.y);
  image.Context()->stroke();

  p = GetPoint(image, n + 2, z, false);
  image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5, "Closed neutron shell",
      CairoImage::HorizontalTextPosition::Left);

  z -= dz;
  p0 = GetPoint(image, n, z, true);
  p0.x -= Opt::lineWidth / 2.0;
  p1 = GetPoint(image, n + 1, z + 1, true);

  image.Context()->move_to(p0.x, p0.y);
  image.Context()->line_to(p1.x, p0.y);
  image.Context()->stroke();

  image.Context()->move_to(p1.x, p1.y);
  image.Context()->line_to(p0.x, p1.y);
  image.Context()->stroke();

  p = GetPoint(image, n + 2, z, false);
  image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5, "Closed proton shell",
      CairoImage::HorizontalTextPosition::Left);

  // make color bar title
  font = Cairo::ToyFontFace::create(Opt::fontFamily, Cairo::FONT_SLANT_NORMAL,
      Cairo::FONT_WEIGHT_BOLD);
  image.Context()->set_font_face(font);

  image.DrawText(pLowerLeft.x - width / 2.0, pUpperRight.y - Opt::fontSize,
      "Abundance", CairoImage::HorizontalTextPosition::Center,
      CairoImage::VerticalTextPosition::Bottom);

  /*
  // write mass fraction
  image.DrawText(pLowerLeft.x,
      pUpperRight.y - width - Opt::fontSize,
      "Mass", CairoImage::HorizontalTextPosition::Center,
      CairoImage::VerticalTextPosition::Bottom);
  image.DrawText(pLowerLeft.x,
      pUpperRight.y - width + 0.2 * Opt::fontSize,
      "fraction", CairoImage::HorizontalTextPosition::Center,
      CairoImage::VerticalTextPosition::Bottom);
  */

  // make arrow
  double length = 75.0;
  double x = 75.0;
  double y = 0.0;

  image.DrawArrow(x, image.Height() - 1.0 - y,
      x + length, image.Height() - 1.0 - y);
  image.DrawText(x + length + Opt::fontSize / 2.0, image.Height() - 1.0 - y, "N",
      CairoImage::HorizontalTextPosition::Left);

  x = 0.0;
  y = 60.0;
  image.DrawArrow(x, image.Height() - 1.0 - y,
      x, image.Height() - 1.0 - y - length);
  image.DrawText(x, image.Height() - 1.0 - y - length - Opt::fontSize / 2.0,
      "Z", CairoImage::HorizontalTextPosition::Center,
      CairoImage::VerticalTextPosition::Bottom);


  font = Cairo::ToyFontFace::create(
        Opt::fontFamily, Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    image.Context()->set_font_face(font);
    image.Context()->set_font_size(Opt::fontSize);
}

void NucleiChart::DrawNumbers(CairoImage& image) const {
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  Cairo::RefPtr<Cairo::ToyFontFace> font = Cairo::ToyFontFace::create(
      Opt::fontFamily, Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
  image.Context()->set_font_face(font);
  image.Context()->set_font_size(10.0);

  for (int n = 0; n <= mDat.MaxN(); n += 4) {
    int offset = 2;
    int top = mDat.ZLimitsVsN()[n].min - offset;

    // move the number furhter down if it's too close to the border from the
    // columns to the left or right
    if (n > 0) {
      top = std::min(top, mDat.ZLimitsVsN()[n - 1].min - offset / 2);
      if (n > 1)
        top = std::min(top, mDat.ZLimitsVsN()[n - 2].min - offset / 2);
    }

    if (n < mDat.MaxN()) {
      top = std::min(top, mDat.ZLimitsVsN()[n + 1].min - offset / 2);
      if (n < mDat.MaxN() - 1)
        top = std::min(top, mDat.ZLimitsVsN()[n + 2].min - offset / 2);
    }

    Point p = GetPoint(image, n, top, false);
    // move half the cell size to the right and half a pixel to the left,
    // so that the center is at the center of the cell
    image.DrawText(p.x + Opt::cellSize / 2.0 - 0.5, p.y,
        std::to_string(n), CairoImage::HorizontalTextPosition::Center,
        CairoImage::VerticalTextPosition::Top);
  }

  for (int z = 0; z <= mDat.MaxZ(); z += 2) {
    int offset = 2;
    int left = mDat.NLimitsVsZ()[z].min - offset;

    // move the number further to the left if it's too close to the border from
    // the row below or above
    if (z > 0) {
      left = std::min(left, mDat.NLimitsVsZ()[z - 1].min - offset / 2);
    }
    if (z < mDat.MaxZ()) {
      left = std::min(left, mDat.NLimitsVsZ()[z + 1].min - offset / 2);
    }

    Point p = GetPoint(image, left, z, false);
    // move half the cell size to the right and half a pixel to the left,
    // so that the center is at the center of the cell
    image.DrawText(p.x, p.y - Opt::cellSize / 2.0 - 0.5,
        std::to_string(z), CairoImage::HorizontalTextPosition::Right,
        CairoImage::VerticalTextPosition::Middle);
  }
}

void NucleiChart::DrawMissing(CairoImage& image) const {
  // TODO make this automatic
  DrawNuclide(image, 4, 5, Opt::missingFillRGBA);
  DrawNuclide(image, 4, 4, Opt::missingFillRGBA);
  DrawNuclide(image, 185, 90, Opt::missingFillRGBA);
  DrawNuclide(image, 186, 90, Opt::missingFillRGBA);
  DrawNuclide(image, 187, 90, Opt::missingFillRGBA);
  DrawNuclide(image, 188, 90, Opt::missingFillRGBA);
  DrawNuclide(image, 220, 108, Opt::missingFillRGBA);
  DrawNuclide(image, 219, 109, Opt::missingFillRGBA);
  DrawNuclide(image, 225, 108, Opt::missingFillRGBA);
  DrawNuclide(image, 226, 107, Opt::missingFillRGBA);
  DrawNuclide(image, 226, 106, Opt::missingFillRGBA);
  DrawNuclide(image, 226, 105, Opt::missingFillRGBA);
  DrawNuclide(image, 226, 104, Opt::missingFillRGBA);
  DrawNuclide(image, 226, 103, Opt::missingFillRGBA);
  DrawNuclide(image, 227, 105, Opt::missingFillRGBA);
  DrawNuclide(image, 227, 104, Opt::missingFillRGBA);
  DrawNuclide(image, 227, 103, Opt::missingFillRGBA);
}

Point NucleiChart::GetPoint(CairoImage& image, const int n, const int z,
    const bool borderCentered) const {
  double x = (Opt::cellSize + Opt::lineWidth / 2.0) * double(n) +
      Opt::lineWidth / 2.0;
  double y = (Opt::cellSize + Opt::lineWidth / 2.0) * double(z);

  if (borderCentered) {
    x -= Opt::lineWidth / 4.0;
    y -= Opt::lineWidth / 4.0;
  }

  return Point(x, image.Height() - 1.0 - y);
}

Point NucleiChart::GetCellCenterPoint(CairoImage& image, const int n, const int z,
    const bool borderCentered) const {
  auto lowerLeft = GetPoint(image, n, z, borderCentered);
  auto upperRight = GetPoint(image, n + 1, z + 1, borderCentered);

  return Point(0.5 * (lowerLeft.x + upperRight.x),
      0.5 * (lowerLeft.y + upperRight.y));
}


void NucleiChart::DrawMagicNumberLine(CairoImage& image, const int constNum,
    const int low, const int high, const bool constN) const {
  Point p0(0, 0);
  Point p1(0, 0);
  if (constN) {
    p0 = GetPoint(image, constNum, low, true);
    p1 = GetPoint(image, constNum + 1, high, true);
  } else {
    p0 = GetPoint(image, low, constNum, true);
    p1 = GetPoint(image, high, constNum + 1, true);
  }

  auto rgb = Opt::magicNumberLineRGB;
  image.Context()->set_source_rgb(rgb[0], rgb[1], rgb[2]);

  image.Context()->move_to(p0.x, p0.y);
  if (constN)
    image.Context()->line_to(p0.x, p1.y);
  else
    image.Context()->line_to(p1.x, p0.y);
  image.Context()->stroke();

  image.Context()->move_to(p1.x, p1.y);
  if (constN)
    image.Context()->line_to(p1.x, p0.y);
  else
    image.Context()->line_to(p0.x, p1.y);
  image.Context()->stroke();
}

void NucleiChart::DrawMagicNumberAndPeakLines(CairoImage& image) const {
  image.Context()->set_line_width(0.75 * Opt::lineWidth);
  for (int n : Opt::magicN) {
    DrawMagicNumberLine(image, n,
        mDat.ZLimitsVsN()[n].min - Opt::magicNLineOvershoot,
        mDat.ZLimitsVsN()[n].max + 1 + Opt::magicNLineOvershoot, true);
  }

  for (int z : Opt::magicN) {
    DrawMagicNumberLine(image, z,
        mDat.NLimitsVsZ()[z].min - Opt::magicNLineOvershoot,
        mDat.NLimitsVsZ()[z].max + 1 + Opt::magicNLineOvershoot, false);
  }

  auto rgb = Opt::peakLineRGB;
  image.Context()->set_source_rgb(rgb[0], rgb[1], rgb[2]);

  for (int a : Opt::leftBinEdges) {
    if (a < 10)
      continue;
    DrawConstALine(image, a);
  }

  image.Context()->set_line_width(Opt::lineWidth);
}

void NucleiChart::DrawConstALine(CairoImage& image, const int a) const {

  int nMin = mDat.NLimitsVsA()[a].min - Opt::peakLineOvershoot;
  int nMax = mDat.NLimitsVsA()[a].max + Opt::peakLineOvershoot;

  Point p0 = GetPoint(image, nMin, a - nMin + 1, true);
  image.Context()->move_to(p0.x, p0.y);

  for (int n = nMin; n <= nMax; ++n) {
    int z = a - n;
    Point p0 = GetPoint(image, n, z + 1, true);
    Point p1 = GetPoint(image, n + 1, z + 1, true);
    Point p2 = GetPoint(image, n + 1, z, true);

    image.Context()->line_to(p0.x, p0.y);
    image.Context()->line_to(p1.x, p1.y);
    image.Context()->line_to(p2.x, p2.y);
  }

  image.Context()->stroke();
}

void NucleiChart::DrawPathLines(CairoImage& image, const int timeIdx) const {
  if (!mDat.HavePathInfo())
    return;

  // find max Y/tau_beta for each Z
  std::vector<double> maxYDivTauBeta(mDat.MaxZ() + 1, 1.0e-50);
  std::vector<int> NOfMaxYDivTauBeta(mDat.MaxZ() + 1, -1);

  for (int i = 0; i < mDat.NumIsotopes(); ++i) {
    int z = mDat.Z()[i];
    double yDivTauBeta = mDat.YDivTauBetaVsIsotopeVsTLogInterp()[i][timeIdx];

    if (yDivTauBeta > maxYDivTauBeta[z]) {
      maxYDivTauBeta[z] = yDivTauBeta;
      NOfMaxYDivTauBeta[z] = mDat.N()[i];
    }
  }

  image.Context()->save();
  image.Context()->set_line_width(Opt::lineWidth);
  image.Context()->set_antialias(Cairo::ANTIALIAS_NONE);
  auto rgb = Opt::pathLineRGB;
  image.Context()->set_source_rgb(rgb[0], rgb[1], rgb[2]);

  // find first z with a maxYDivTauBeta
  for (int z = 0; z < mDat.MaxZ(); ++z) {
    if ((NOfMaxYDivTauBeta[z] > -1) && (NOfMaxYDivTauBeta[z+1] > -1)) {
      Point p0 = GetCellCenterPoint(image, NOfMaxYDivTauBeta[z], z, true);
      Point p1 = GetCellCenterPoint(image, NOfMaxYDivTauBeta[z] - 1, z+1, true);
      Point p2 = GetCellCenterPoint(image, NOfMaxYDivTauBeta[z+1], z+1, true);

      image.Context()->move_to(p0.x, p0.y);
      image.Context()->line_to(p1.x, p1.y);
      image.Context()->line_to(p2.x, p2.y);
      image.Context()->stroke();
    }
  }

  image.Context()->restore();
}

void NucleiChart::DrawStableNuclide(CairoImage& image, const int N,
    const int Z) const {
  DrawCell(image, N, Z, {{ 0.0, 0.0, 0.0, 0.0 }}, Opt::stableBorderRGB, true);
}

void NucleiChart::DrawStableNuclides(CairoImage& image) const {
  for (unsigned int i = 0; i < Opt::stableN.size(); ++i) {
    DrawStableNuclide(image, Opt::stableN[i], Opt::stableZ[i]);
  }
}

void NucleiChart::DrawCell(CairoImage& image, const int N, const int Z,
    const std::array<double, 4> fillRGBA,
    const std::array<double, 3> borderRGB, const bool drawBorderOnly) const {

  Point pt = GetPoint(image, N, Z, false);
  double size = Opt::cellSize;

  // fill
  if (!drawBorderOnly) {
    image.Context()->set_source_rgba(fillRGBA[0], fillRGBA[1], fillRGBA[2],
        fillRGBA[3]);
    image.Context()->rectangle(pt.x, pt.y, size, -size);
    image.Context()->fill();
  }

  // border
  pt.y += Opt::lineWidth / 2.0;
  size += Opt::lineWidth / 2.0;
  image.Context()->set_source_rgb(borderRGB[0], borderRGB[1], borderRGB[2]);
  image.Context()->rectangle(pt.x, pt.y, size, -size);
  image.Context()->stroke();
}

std::vector<double> NucleiChart::DrawNuclides(CairoImage& image,
    const int timeIdx, const bool drawBackground) const {
  image.Context()->save();
  image.Context()->set_line_width(Opt::lineWidth / 2.0);
  image.Context()->set_antialias(Cairo::ANTIALIAS_NONE);

  std::vector<double> binCounts(Opt::leftBinEdges.size());

  for (int i = 0; i < mDat.NumIsotopes(); ++i) {
    int a = mDat.A()[i];
    double massAbund = mDat.MassAbundanceVsIsotopeVsTLinInterp()[i][timeIdx];
    double abund = mDat.AbundanceVsIsotopeVsTLogInterp()[i][timeIdx];

    if (!drawBackground && (abund == 0.0))
      continue;

    // find bin
    int bin = Opt::leftBinEdges.size() - 1;
    for (unsigned int j = 1; j < Opt::leftBinEdges.size(); ++j) {
      if (a <= Opt::leftBinEdges[j]) {
        bin = j - 1;
        break;
      }
    }
    binCounts[bin] += massAbund;

    if (drawBackground)
      DrawNuclide(image, mDat.N()[i], mDat.Z()[i], Opt::binColors[bin]);
    else
      DrawNuclide(image, mDat.N()[i], mDat.Z()[i],
          mpColorScale->GetRGBA(abund));
  }

  DrawMissing(image);
  DrawStableNuclides(image);
  image.Context()->restore();

  return binCounts;
}

void NucleiChart::DrawPlotDots(CairoImage& image, const int timeIdx) const {
  if (mPlotCoords.size() == 0)
    return;

  image.Context()->save();
  image.Context()->set_identity_matrix();
  image.Context()->translate(mXOffset / 2.0, mHeight - mChartHeight + mYOffset
      - 5);

  int idx = timeIdx - mDat.FirstTimeIndex();

  image.Context()->arc(mPlotCoords[idx].t,
      mPlotHeight - mPlotCoords[idx].temp, 4.0, 0.0, 2.0 * M_PI);
  image.Context()->set_source_rgb(1.0, 0.0, 0.0);
  image.Context()->fill();


  image.Context()->arc(mPlotCoords[idx].t,
      mPlotHeight - mPlotCoords[idx].density, 4.0, 0.0, 2.0 * M_PI);
  image.Context()->set_source_rgb(0.0, 0.0, 1.0);
  image.Context()->fill();

  image.Context()->arc(mPlotCoords[idx].t,
      mPlotHeight - mPlotCoords[idx].edot, 4.0, 0.0, 2.0 * M_PI);
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  image.Context()->fill();

  image.Context()->restore();
}

void NucleiChart::DrawText(CairoImage& image, const int timeIdx) const {
  char buf[128];
  image.Context()->set_source_rgb(0.0, 0.0, 0.0);
  Cairo::RefPtr<Cairo::ToyFontFace> font = Cairo::ToyFontFace::create(
      Opt::fontFamily, Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
  image.Context()->set_font_face(font);
  image.Context()->set_font_size(Opt::fontSize);

  // TODO
  double x = 460.0;
  double y = image.Height() - mChartHeight;
  double ySkip = 28.0;

  sprintf(buf, "%1.2E K", mDat.TemperatureVsT()[timeIdx] * 1.0E9);
  image.DrawText(x, y, "Temperature = " +
      CairoImage::ReplaceMinusWithEndash(std::string(buf)),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Top);

  y += ySkip;
  sprintf(buf, "%1.2E g / cm^{3}", mDat.RhoVsT()[timeIdx]);
  image.DrawText(x, y, "Density = " +
      CairoImage::ReplaceMinusWithEndash(std::string(buf)),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Top);

  y += ySkip;
  sprintf(buf, "%1.2E erg / s / g",
      /*mrDat.Mass() * */mDat.EdotVsT()[timeIdx]);
  image.DrawText(x, y, "Heating rate = " +
      CairoImage::ReplaceMinusWithEndash(std::string(buf)),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Top);

  y += ySkip;
  sprintf(buf, "%1.2E kB / baryon", mDat.EntropyVsT()[timeIdx]);
  image.DrawText(x, y, "Entropy = " +
      CairoImage::ReplaceMinusWithEndash(std::string(buf)),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Top);

  y += ySkip;
  sprintf(buf, "%1.3f", mDat.YeVsT()[timeIdx]);
  image.DrawText(x, y, "Ye = " +
      CairoImage::ReplaceMinusWithEndash(std::string(buf)),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Top);


  // time

  int n = 136;
  int z = 40;
  Point p = GetPoint(image, n, z, false);

  x = p.x;
  y = p.y - Opt::cellSize / 2.0 - 0.5;
  ySkip = 25.0;

  font = Cairo::ToyFontFace::create(
      Opt::fontFamily, Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_BOLD);
  image.Context()->set_font_face(font);

  double t = mDat.Time()[timeIdx];
  sprintf(buf, "%1.2E s", t);
  std::string timeInS(buf);

  double humanT = t;
  std::string unit = "?";

  if (t >= 1.0E-9) {
    humanT = t * 1.0E9;
    unit = "ns";
  }
  if (t >= 1.0E-6) {
    humanT = t * 1.0E6;
    unit = "\u00B5s";
  }
  if (t >= 1.0E-3) {
    humanT = t * 1.0E3;
    unit = "ms";
  }
  if (t >= 1.0) {
    humanT = t;
    unit = "s";
  }
  if (t >= 60.0) {
    humanT = t / 60.0;
    unit = "min";
  }
  if (t >= (60.0 * 60.0)) {
    humanT = t / (60.0 * 60.0);
    unit = "hr";
  }
  if (t >= (24.0 * 60.0 * 60.0)) {
    humanT = t / (24.0 * 60.0 * 60.0);
    unit = "day";
  }
  if (t >= (365.0 * 24.0 * 60.0 * 60.0)) {
    humanT = t / (365.0 * 24.0 * 60.0 * 60.0);
    unit = "yr";
  }

  if (humanT >= 100.0)
    humanT = 10.0 * round(humanT / 10.0);

  if (humanT < 10.0)
    sprintf(buf, "%.1f", humanT);
  else
    sprintf(buf, "%.0f", humanT);

  std::string humanTStr(buf);

  Cairo::TextExtents extents;
  image.Context()->get_text_extents("Time", extents);

  image.DrawText(x, y, "Time",
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Bottom);
  image.DrawText(x + extents.width + 3.0, y, " = " +
      CairoImage::ReplaceMinusWithEndash(timeInS),
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Bottom);

  y += ySkip;
  image.DrawText(x + extents.width + 3.0, y, " = " + humanTStr + " " + unit,
      CairoImage::HorizontalTextPosition::Left,
      CairoImage::VerticalTextPosition::Bottom);
}

