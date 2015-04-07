/// \file Scale.hpp
/// \author jonas
/// \since Jun 14, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_SCALE_HPP_
#define SKYNET_MOVIEMAKER_SCALE_HPP_

#ifdef SKYNET_USE_MOVIE

#include <string>
#include <vector>

class Scale {
public:
  Scale(const double xMin, const double xMax, const bool useLog,
      const std::vector<double>& tics,
      const std::vector<std::string>& ticLabels);

  double operator()(const double x) const;

  double XMin() const {
    return mXMin;
  }

  double XMax() const {
    return mXMax;
  }

  double Dx() const {
    return mDx;
  }

  bool UseLog() const {
    return mUseLog;
  }

  const std::vector<double>& Tics() const {
    return mTics;
  }

  const std::vector<std::string>& TicLabels() const {
    return mTicLabels;
  }

private:
  double mXMin;
  double mXMax;
  double mDx;
  bool mUseLog;

  std::vector<double> mTics;
  std::vector<std::string> mTicLabels;
};

#endif // SKYNET_USE_MOVIE

#endif // SKYNET_MOVIEMAKER_SCALE_HPP_
