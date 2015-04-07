/// \file PiecewiseLinearColorScale.hpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_PIECEWISELINEARCOLORSCALE_HPP_
#define SKYNET_MOVIEMAKER_PIECEWISELINEARCOLORSCALE_HPP_

#ifdef SKYNET_USE_MOVIE

#include <array>
#include <vector>

#include "MovieMaker/ColorScale.hpp"
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class PieceWiseLinearColorScale : public ColorScale {
public:
  PieceWiseLinearColorScale(const double xMin, const double xMax,
      const bool useLog, const std::vector<double>& tics,
      const std::vector<std::string>& ticLabels,
      const std::vector<std::pair<double, double>>& red,
      const std::vector<std::pair<double, double>>& green,
      const std::vector<std::pair<double, double>>& blue,
      const std::vector<std::pair<double, double>>& alpha);

  virtual std::array<double, 4> GetRGBARange01(const double x) const;

private:
  PiecewiseLinearFunction mRed;
  PiecewiseLinearFunction mGreen;
  PiecewiseLinearFunction mBlue;
  PiecewiseLinearFunction mAlpha;
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_PIECEWISELINEARCOLORSCALE_HPP_
