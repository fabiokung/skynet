/// \file ColorScale.hpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_COLORSCALE_HPP_
#define SKYNET_MOVIEMAKER_COLORSCALE_HPP_

#ifdef SKYNET_USE_MOVIE

#include <array>
#include <string>
#include <vector>

#include "MovieMaker/Scale.hpp"

class ColorScale {
public:
  ColorScale(const double xMin, const double xMax, const bool useLog,
      const std::vector<double>& tics,
      const std::vector<std::string>& ticLabels);

  virtual ~ColorScale() { }

  std::array<double, 4> GetRGBA(const double x) const;

  virtual std::array<double, 4> GetRGBARange01(const double x) const =0;

  const Scale& GetScale() const {
    return mScale;
  }

private:
  Scale mScale;
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_COLORSCALE_HPP_
