/// \file CairoImage.hpp
/// \author jlippuner
/// \since Jun 12, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_CAIROIMAGE_HPP_
#define SKYNET_MOVIEMAKER_CAIROIMAGE_HPP_

#include <array>
#include <string>

#include <cairomm/context.h>
#include <cairomm/surface.h>


class CairoImage {
public:
  enum class HorizontalTextPosition {
    Left,
    Center,
    Right
  };

  enum class VerticalTextPosition {
    Top,
    Middle,
    Bottom
  };

  CairoImage(const int width, const int height,
      const std::array<double, 3>& backgroundRGB);

  CairoImage(const int width, const int height):
    CairoImage(width, height, {{0.0, 0.0, 0.0}}) { }

  void Save(const std::string& path);

  void SetPixel(const int x, const int y, const std::array<double, 3>& rgb);

  static std::string ReplaceMinusWithEndash(const std::string& str);

  void DrawText(const double x, const double y,
      const std::vector<std::string>& strs, const std::vector<bool>& super,
      const HorizontalTextPosition hPos = HorizontalTextPosition::Center,
      const VerticalTextPosition vPos = VerticalTextPosition::Middle);

  void DrawText(const double x, const double y, const std::string& text,
      const HorizontalTextPosition hPos = HorizontalTextPosition::Center,
      const VerticalTextPosition vPos = VerticalTextPosition::Middle);

  void DrawArrow(const double start_x, const double start_y,
      const double end_x, const double end_y) const;

  const Cairo::RefPtr<Cairo::ImageSurface>& Surface() const {
    return mSurface;
  }

  Cairo::RefPtr<Cairo::ImageSurface>& Surface() {
    return mSurface;
  }

  const Cairo::RefPtr<Cairo::Context>& Context() const {
    return mContext;
  }

  Cairo::RefPtr<Cairo::Context>& Context() {
    return mContext;
  }

  double Height() const {
    return double(mSurface->get_height());
  }

private:
  Cairo::RefPtr<Cairo::ImageSurface> mSurface;
  Cairo::RefPtr<Cairo::Context> mContext;
};

#endif // SKYNET_MOVIEMAKER_CAIROIMAGE_HPP_
