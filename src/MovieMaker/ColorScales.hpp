/// \file ColorScales.hpp
/// \author jlippuner
/// \since Jun 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_COLORSCALES_HPP_
#define SKYNET_MOVIEMAKER_COLORSCALES_HPP_

#ifdef SKYNET_USE_MOVIE

#include "MovieMaker/PiecewiseLinearColorScale.hpp"

class WhiteBlueGreenYellowRedColorScale: public PieceWiseLinearColorScale {
public:
  WhiteBlueGreenYellowRedColorScale(const double xMin, const double xMax,
      const bool useLog, const std::vector<double>& tics,
      const std::vector<std::string>& ticLabels) :
      PieceWiseLinearColorScale(xMin, xMax, useLog, tics, ticLabels,
          {
            {   0.0 / 253.0, 255.0 / 255.0 },
            {  36.0 / 253.0, 157.0 / 255.0 },
            {  72.0 / 253.0,  72.0 / 255.0 },
            { 108.0 / 253.0,  73.0 / 255.0 },
            { 145.0 / 253.0, 250.0 / 255.0 },
            { 181.0 / 253.0, 245.0 / 255.0 },
            { 217.0 / 253.0, 211.0 / 255.0 },
            { 253.0 / 253.0, 146.0 / 255.0 }
          },
          {
            {   0.0 / 253.0, 255.0 / 255.0 },
            {  36.0 / 253.0, 218.0 / 255.0 },
            {  72.0 / 253.0, 142.0 / 255.0 },
            { 108.0 / 253.0, 181.0 / 255.0 },
            { 145.0 / 253.0, 232.0 / 255.0 },
            { 181.0 / 253.0, 106.0 / 255.0 },
            { 217.0 / 253.0,  31.0 / 255.0 },
            { 253.0 / 253.0,  21.0 / 255.0 }
          },
          {
            {   0.0 / 253.0, 255.0 / 255.0 },
            {  36.0 / 253.0, 247.0 / 255.0 },
            {  72.0 / 253.0, 202.0 / 255.0 },
            { 108.0 / 253.0,  70.0 / 255.0 },
            { 145.0 / 253.0,  92.0 / 255.0 },
            { 181.0 / 253.0,  41.0 / 255.0 },
            { 217.0 / 253.0,  40.0 / 255.0 },
            { 253.0 / 253.0,  25.0 / 255.0 }
          },
          {
            { 0.0, 0.0 },
            { 18.0 / 253.0, 1.0 },
            { 1.0, 1.0 }
          }) { }
};


class WhiteBlackColorScale: public PieceWiseLinearColorScale {
public:
  WhiteBlackColorScale(const double xMin, const double xMax,
      const bool useLog, const std::vector<double>& tics,
      const std::vector<std::string>& ticLabels) :
      PieceWiseLinearColorScale(xMin, xMax, useLog, tics, ticLabels,
          {
            { 0.0, 1.0 },
            { 1.0, 0.0 }
          },
          {
            { 0.0, 1.0 },
            { 1.0, 0.0 }
          },
          {
            { 0.0, 1.0 },
            { 1.0, 0.0 }
          },
          {
            { 0.0, 0.0 },
            { 0.1, 1.0 },
            { 1.0, 1.0 }
          }) { }
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_COLORSCALES_HPP_
