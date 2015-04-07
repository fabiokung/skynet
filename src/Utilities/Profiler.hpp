/// \file Profiler.hpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_PROFILER_HPP_
#define SKYNET_UTILITIES_PROFILER_HPP_

#include <string>

#include "Network/NetworkOutput.hpp"
#include "Utilities/StopWatch.hpp"

class Profiler {
public:
  StopWatch Everything;
  StopWatch RateCalculation;
  StopWatch YdotCalculation;
  StopWatch JacobianCalculation;
  StopWatch MatrixInversion;
  StopWatch NSEEvolution;
//  StopWatch MatrixInversionCuda;
//  StopWatch MatrixInversionCudaDataCopying;
//  StopWatch MatrixInversionCudaWithoutDataCopying;
//  StopWatch Cusparse;

  void Report(NetworkOutput * const pOutput, const int numNRIterations);

private:
  void Report(NetworkOutput * const pOutput, const std::string& description,
      const StopWatch& stopWatch, const bool subtractFromTotal = true);

  void Report(NetworkOutput * const pOutput, const std::string& description,
      const double elapsedTimeInSeconds, const bool subtractFromTotal = true);

  double mTotal;
  double mOther;

  double mNumTimeSteps;
  double mNumNRIterations;
};

#endif // SKYNET_UTILITIES_PROFILER_HPP_
