/// \file Profiler.cpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#include "Utilities/Profiler.hpp"

#include <stdexcept>
#include <stdio.h> 

void Profiler::Report(NetworkOutput * const pOutput, const int numNRIterations) {
  if (Everything.IsRunning())
    throw std::logic_error("Stop watch Everything still running when "
        "reporting.");
  if (RateCalculation.IsRunning())
    throw std::logic_error("Stop watch RateCalculation still running when "
        "reporting.");
  if (YdotCalculation.IsRunning())
    throw std::logic_error("Stop watch YdotCalculation still running when "
        "reporting.");
  if (JacobianCalculation.IsRunning())
    throw std::logic_error("Stop watch JacobianCalculation still running when "
        "reporting.");
  if (MatrixInversion.IsRunning())
    throw std::logic_error("Stop watch MatrixInversion still running when "
        "reporting.");
  if (NSEEvolution.IsRunning())
    throw std::logic_error("Stop watch NSEEvolution still running when "
        "reporting.");
//  if (MatrixInversionCuda.IsRunning())
//    throw std::logic_error("Stop watch MatrixInversionCuda still running when "
//        "reporting.");
//  if (MatrixInversionCudaDataCopying.IsRunning())
//    throw std::logic_error("Stop watch MatrixInversionCudaDataCopying still running when "
//        "reporting.");
//  if (MatrixInversionCudaWithoutDataCopying.IsRunning())
//    throw std::logic_error("Stop watch MatrixInversionCudaWithoutDataCopying still running when "
//        "reporting.");
//  if (Cusparse.IsRunning())
//    throw std::logic_error("Stop watch Cusparse still running when "
//        "reporting.");

  mTotal = Everything.GetElapsedTimeInSeconds();
  mOther = mTotal;

  mNumTimeSteps = (double)pOutput->NumEntries();
  mNumNRIterations = (double)numNRIterations;

  pOutput->Log("# Profiling report:\n");
  pOutput->Log("#\n");
  pOutput->Log("# %lu time steps, %i Newton-Raphson iterations\n",
      pOutput->NumEntries(), numNRIterations);
  pOutput->Log("#\n");
  pOutput->Log("# %16s  %12s  %13s  %13s  %8s\n", "", "total", "per time step",
      "per NR iter", "fraction");
  Report(pOutput, "calc rates", RateCalculation);
  Report(pOutput, "calc Ydot", YdotCalculation);
  Report(pOutput, "calc Jacobian", JacobianCalculation);
  Report(pOutput, "matrix inversion", MatrixInversion);
  Report(pOutput, "NSE evolution", NSEEvolution);
//  Report(pOutput, "matrix inv CUDA", MatrixInversionCuda);
//  Report(pOutput, "mat inv CUDA copying", MatrixInversionCudaDataCopying, false);
//  Report(pOutput, "mat inv CUDA w/o copy", MatrixInversionCudaWithoutDataCopying, false);
//  Report(pOutput, "cusparse", Cusparse, false);
  Report(pOutput, "other", mOther);
  Report(pOutput, "total", mTotal);
}

void Profiler::Report(NetworkOutput * const pOutput,
    const std::string& description, const StopWatch& stopWatch,
    const bool subtractFromTotal) {
  Report(pOutput, description, stopWatch.GetElapsedTimeInSeconds(),
      subtractFromTotal);
}

void Profiler::Report(NetworkOutput * const pOutput,
    const std::string& description, const double elapsedTimeInSeconds,
    const bool subtractFromTotal) {
  if (subtractFromTotal)
    mOther -= elapsedTimeInSeconds;

  pOutput->Log("# %16s", description.c_str());
  pOutput->Log("  %10.2f s  %10.2f ms  %10.2f ms  %7.2f%%\n",
      elapsedTimeInSeconds, 1000.0 * elapsedTimeInSeconds / (mNumTimeSteps + 1.e-10),
      1000.0 * elapsedTimeInSeconds / (mNumNRIterations + 1.e-10),
      100.0 * elapsedTimeInSeconds / (mTotal + 1.e-10));
}
