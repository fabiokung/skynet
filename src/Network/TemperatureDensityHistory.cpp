/// \file TemperatureDensityHistory.cpp
/// \author jlippuner
/// \since Aug 11, 2014
///
/// \brief
///
///

#include <Utilities/File.hpp>
#include "Network/TemperatureDensityHistory.hpp"

#include <fstream>


TemperatureDensityHistory TemperatureDensityHistory::CreateFromValues(
    const double Ye, const double startTime, const double endTime,
    const std::vector<double>& times, const std::vector<double>& temperatures,
    const std::vector<double>& densities, const bool interpolateInLogSpace) {
  return TemperatureDensityHistory(Ye, startTime, endTime, times, temperatures,
      densities, interpolateInLogSpace);
}

TemperatureDensityHistory TemperatureDensityHistory::CreateFromFile(
    const std::string& filePath, const bool interpolateInLogSpace) {
  if (!File::Exists(filePath)) {
    throw std::invalid_argument("The temperature and density history file '"
        + filePath + "' does not exist.");
  }

  double Ye, t0, tFinal;
  std::vector<double> ts, T9s, rhos;

  std::ifstream ifs(filePath.c_str(), std::ifstream::in);
  ifs >> Ye >> t0 >> tFinal;

  while (!ifs.eof()) {
    double t, T9, rho;
    ifs >> t >> T9 >> rho;

    if (!ifs.good())
      break;

    ts.push_back(t);
    T9s.push_back(T9);
    rhos.push_back(rho);
  }

  ifs.close();

  return CreateFromValues(Ye, t0, tFinal, ts, T9s, rhos, interpolateInLogSpace);
}

TemperatureDensityHistory::TemperatureDensityHistory(const double Ye,
    const double startTime, const double endTime,
    const std::vector<double>& times, const std::vector<double>& temperatures,
    const std::vector<double>& densities, const bool interpolateInLogSpace) :
    mYe(Ye),
    mStartTime(startTime),
    mEndTime(endTime),
    mTimes(times),
    mTemperatures(temperatures),
    mDensities(densities),
    mTemperatureVsTime(mTimes, mTemperatures, interpolateInLogSpace),
    mDensityVsTime(mTimes, mDensities, interpolateInLogSpace) { }
