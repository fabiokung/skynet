/// \file TemperatureDensityHistory.hpp
/// \author jlippuner
/// \since Aug 11, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_TEMPERATUREDENSITYHISTORY_HPP_
#define SKYNET_NETWORK_TEMPERATUREDENSITYHISTORY_HPP_

#include <string>
#include <vector>

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class TemperatureDensityHistory {
public:
  static TemperatureDensityHistory CreateFromValues(const double Ye,
      const double startTime, const double endTime,
      const std::vector<double>& times, const std::vector<double>& temperatures,
      const std::vector<double>& densities,
      const bool interpolateInLogSpace = true);

  static TemperatureDensityHistory CreateFromFile(const std::string& filePath,
      const bool interpolateInLogSpace = true);

  const std::vector<double>& Densities() const {
    return mDensities;
  }

  double EndTime() const {
    return mEndTime;
  }

  double StartTime() const {
    return mStartTime;
  }

  const std::vector<double>& Temperatures() const {
    return mTemperatures;
  }

  const std::vector<double>& Times() const {
    return mTimes;
  }

  double Ye() const {
    return mYe;
  }

  const PiecewiseLinearFunction& DensityVsTime() const {
    return mDensityVsTime;
  }

  const PiecewiseLinearFunction& TemperatureVsTime() const {
    return mTemperatureVsTime;
  }

private:
  TemperatureDensityHistory(const double Ye, const double startTime,
      const double endTime, const std::vector<double>& times,
      const std::vector<double>& temperatures,
      const std::vector<double>& densities,
      const bool interpolateInLogSpace);

  // raw data
  const double mYe;
  const double mStartTime; // in second
  const double mEndTime; // in second
  const std::vector<double> mTimes; // in second
  const std::vector<double> mTemperatures; // in Giga Kelvin
  const std::vector<double> mDensities; // in gram / cm^3

  // temperature and density vs time functions (linearly interpolated)
  const PiecewiseLinearFunction mTemperatureVsTime;
  const PiecewiseLinearFunction mDensityVsTime;
};

#endif // SKYNET_NETWORK_TEMPERATUREDENSITYHISTORY_HPP_
