/// \file CairoImage.cpp
/// \author jlippuner
/// \since Jun 12, 2014
///
/// \brief
///
///

#include "MovieMaker/CairoImage.hpp"

#include <cassert>
#include <stdexcept>

#include "MovieMaker/NucleiChartOptions.hpp"

CairoImage::CairoImage(const int width, const int height,
    const std::array<double, 3>& backgroundRGB):
    mSurface(Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32, width, height)),
    mContext(Cairo::Context::create(mSurface)) {
  mContext->save();
  mContext->set_source_rgb(backgroundRGB[0], backgroundRGB[1],
      backgroundRGB[2]);
  mContext->paint();
  mContext->restore();
}

void CairoImage::Save(const std::string& path) {
  mSurface->write_to_png(path);
}

void CairoImage::SetPixel(const int x, const int y,
    const std::array<double, 3>& rgb) {
  mContext->set_source_rgb(rgb[0], rgb[1], rgb[2]);
  mContext->rectangle(x, mSurface->get_height() - 1 - y, 1, 1);
  mContext->fill();
}

std::string CairoImage::ReplaceMinusWithEndash(const std::string& str) {
  std::string newStr = str;
  std::size_t minusPos = newStr.find("-");
  if (minusPos != std::string::npos)
    newStr.replace(minusPos, 1, u8"\u2013");

  return newStr;
}

void CairoImage::DrawText(const double x, const double y,
    const std::string& text, const HorizontalTextPosition hPos,
    const VerticalTextPosition vPos) {
  std::vector<std::string> strs;
  std::vector<bool> super;

  std::size_t currentPos = 0;
  std::size_t caretPos = text.find("^{", currentPos);
  while (caretPos != std::string::npos) {
    // substring between currentPos and caretPos is standard text
    strs.push_back(text.substr(currentPos, caretPos - currentPos));
    super.push_back(false);

    // substring between { and } is super script
    std::size_t endBracePos = text.find("}", currentPos);
    if ((endBracePos == std::string::npos) || (endBracePos <= caretPos))
      throw std::runtime_error("Unexpected '}' or missing '}'");

    std::string superStr = text.substr(caretPos + 2,
        endBracePos - 2 - caretPos);
    // replace - with endash
    superStr = ReplaceMinusWithEndash(superStr);
    strs.push_back(superStr);
    super.push_back(true);

    currentPos = endBracePos + 1;
    caretPos = text.find("^{", currentPos);
  }

  // the rest of the string is not super script
  strs.push_back(text.substr(currentPos));
  super.push_back(false);

  DrawText(x, y, strs, super, hPos, vPos);
}

void CairoImage::DrawText(const double x, const double y,
    const std::vector<std::string>& strs, const std::vector<bool>& super,
    const HorizontalTextPosition hPos,
    const VerticalTextPosition vPos) {
  // draw dot at text position
  //mContext->arc(x, y, 1.0, 0, 2 * M_PI);
  //mContext->set_source_rgba(0.0, 1.0, 0.0, 0.5);
  //mContext->fill();

  assert(strs.size() == super.size());

  double superScriptMultiplier = 0.65;
  double superScriptOffset = Opt::fontSize * 0.6;

  // get height and width
  double width = 0.0;
  double height = 0.0;

  for (unsigned int i = 0; i < strs.size(); ++i) {
    std::string str = strs[i];
    Cairo::TextExtents extents;

    if (super[i]) {
      mContext->save();
      mContext->set_font_size(Opt::fontSize * superScriptMultiplier);
      mContext->get_text_extents(str, extents);
      //printf("super width of '%s': %f\n", str.c_str(), extents.width);
      //printf("super x_advance of '%s': %f\n", str.c_str(), extents.x_advance);
      mContext->restore();
    }
    else {
      mContext->get_text_extents(str, extents);
      //printf("width of '%s': %f\n", str.c_str(), extents.width);
      //printf("x_advance of '%s': %f\n", str.c_str(), extents.x_advance);
    }

    width += extents.x_advance;

    if (extents.y_bearing < height)
      height = extents.y_bearing;
  }

  //printf("total width: %f\n", width);

  // draw the text
  mContext->move_to(x, y);

  switch (hPos) {
  case HorizontalTextPosition::Left:
    mContext->rel_move_to(0.0, 0.0);
    break;
  case HorizontalTextPosition::Center:
    mContext->rel_move_to(-width / 2.0, 0.0);
    break;
  case HorizontalTextPosition::Right:
    mContext->rel_move_to(-width, 0.0);
    break;
  }

  switch (vPos) {
  case VerticalTextPosition::Bottom:
    mContext->rel_move_to(0.0, 0.0);
    break;
  case VerticalTextPosition::Middle:
    mContext->rel_move_to(0.0, -height / 2.0);
    break;
  case VerticalTextPosition::Top:
    mContext->rel_move_to(0.0, -height);
    break;
  }

  // draw border
  /*
  mContext->save();
  mContext->set_line_width(0.5);
  mContext->rel_line_to(0.0, height);
  mContext->rel_line_to(width, 0.0);
  mContext->rel_line_to(0.0, -height);
  mContext->rel_line_to(-width, 0.0);
  mContext->stroke_preserve();
  mContext->restore();
  */

  for (unsigned int i = 0; i < strs.size(); ++i) {
    std::string str = strs[i];

    if (super[i]) {
      mContext->save();
      mContext->rel_move_to(0.0, -superScriptOffset);
      mContext->set_font_size(Opt::fontSize * superScriptMultiplier);
      mContext->show_text(str);
      mContext->rel_move_to(0.0, superScriptOffset);
      mContext->restore();
    }
    else {
      mContext->show_text(str);
    }
  }
}

void CairoImage::DrawArrow(const double start_x, const double start_y,
    const double end_x, const double end_y) const {
  double arrow_lenght = 12.5;
  double arrow_angle = 20.0 / 180.0 * M_PI;

  double angle = atan2(end_y - start_y, end_x - start_x) + M_PI;

  double x1 = end_x + arrow_lenght * cos(angle - arrow_angle);
  double y1 = end_y + arrow_lenght * sin(angle - arrow_angle);
  double x2 = end_x + arrow_lenght * cos(angle + arrow_angle);
  double y2 = end_y + arrow_lenght * sin(angle + arrow_angle);

  mContext->move_to(start_x, start_y);
  mContext->line_to ((x1 + x2) / 2.0, (y1 + y2) / 2.0);
  mContext->stroke();

  mContext->move_to (x1, y1);
  mContext->line_to (end_x, end_y);
  mContext->line_to (x2, y2);
  mContext->fill();

}
