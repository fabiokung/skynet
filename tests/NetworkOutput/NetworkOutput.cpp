/// \file StopWatch.cpp
/// \author jlippuner
/// \since Jul 10, 2014
///
/// \brief
///
///

#include <cmath>
#include <cstdlib>
#include <vector>

#include "BuildInfo.hpp"
#include "Network/NetworkOptions.hpp"
#include "Network/NetworkOutput.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Utilities/StopWatch.hpp"
#include "Utilities/FloatingPointExceptions.hpp"

struct Entry {
  double time;
  double dt;
  double temperature;
  double density;
  double entropy;
  double heatingRate;
  double Ye;
  std::vector<double> Y;
  std::vector<double> YVsA;
};

Entry MakeEntry(const int timeIdx, const int numNuc, const std::vector<int>& As,
    const int maxA) {
  Entry entry;
  double t = (double)timeIdx;
  entry.time = t;
  entry.dt = -t;
  entry.temperature = 1.3 * t + 13.0;
  entry.density = t * t;
  entry.entropy = sqrt(t);
  entry.heatingRate = 1.0 / (t + 1.0);
  entry.Ye = log(1.0 + t);
  entry.Y = std::vector<double>(numNuc);
  entry.YVsA = std::vector<double>(maxA + 1, 0.0);

  for (int i = 0; i < numNuc; ++i) {
    entry.Y[i] = t * (double)(100 * numNuc) + (double)i;
    entry.YVsA[As[i]] += entry.Y[i];
  }

  return entry;
}

int main(int, char**) {
  FloatingPointExceptions::Enable();

  NuclideLibrary nuclib = NuclideLibrary::CreateFromWebnucleoXML(
      SkyNetRoot + "/data/webnucleo_nuc_v2.0.xml");

  int maxA = 0;
  for (int i = 0; i < nuclib.NumNuclides(); ++i) {
    if (nuclib.As()[i] > maxA)
      maxA = nuclib.As()[i];
  }

  int numNuc = nuclib.NumNuclides();
  int numTimes = 1000;

  // write an output file
  {
    auto output = NetworkOutput::CreateNew("out", nuclib, false);

    std::vector<std::vector<double>> YVsAVsTime;

    for (int t = 0; t < numTimes; ++t) {
      auto e = MakeEntry(t, numNuc, nuclib.As(), maxA);
      YVsAVsTime.push_back(e.YVsA);
      ThermodynamicState thermoState(e.temperature, e.density, e.entropy, 0.0,
          0.0, e.Ye, 0.0);
      output.AddEntry(e.time, { e.dt, "foo" }, 0, 0, thermoState, e.heatingRate,
          e.Y, 1.0);
    }

    // check Y vs A is correct
    for (int t = 0; t < numTimes; ++t) {
      for (int A = 0; A <= maxA; ++A) {
        if (YVsAVsTime[t][A] != output.YVsAVsTime()[t][A]) {
          printf("failed to calculate YVsAVsTime t = %i, A = %i, "
              "%.4E != %.4E\n", t, A, YVsAVsTime[t][A],
              output.YVsAVsTime()[t][A]);
          return EXIT_FAILURE;
        }
      }
    }
  }

  // now read that file from disk and compare to known values
  {
    StopWatch sw;
    sw.Start();
    auto output = NetworkOutput::ReadFromFile("out.h5");
    sw.Stop();
    printf("reading H5 file took %f ms\n", sw.GetElapsedTimeInMilliseconds());

    if (!output.IsReadOnly()) {
      printf("output read from file is not read-only\n");
      return EXIT_FAILURE;
    }

    if (output.NumNuclides() != (std::size_t)numNuc) {
      printf("output read from file has wrong number of nuclides\n");
      return EXIT_FAILURE;
    }

    if (output.NumEntries() != (std::size_t)numTimes) {
      printf("output read from file has wrong number of entries\n");
      return EXIT_FAILURE;
    }

    for (int t = 0; t < numTimes; ++t) {
      auto e = MakeEntry(t, numNuc, nuclib.As(), maxA);

      if (e.time != output.Times()[t])
        throw std::runtime_error("Wrong time at " + std::to_string(t));
      if (e.dt != output.DtVsTime()[t])
        throw std::runtime_error("Wrong dt at " + std::to_string(t));
      if (e.temperature != output.TemperatureVsTime()[t])
        throw std::runtime_error("Wrong temperature at " + std::to_string(t));
      if (e.density != output.DensityVsTime()[t])
        throw std::runtime_error("Wrong density at " + std::to_string(t));
      if (e.entropy != output.EntropyVsTime()[t])
        throw std::runtime_error("Wrong entropy at " + std::to_string(t));
      if (e.heatingRate != output.HeatingRateVsTime()[t])
        throw std::runtime_error("Wrong heating rate at " + std::to_string(t));
      if (e.Ye != output.YeVsTime()[t])
        throw std::runtime_error("Wrong Ye at " + std::to_string(t));

      for (int i = 0; i < numNuc; ++i) {
        if (e.Y[i] != output.YVsTime()[t][i])
          throw std::runtime_error("Wrong Y at " + std::to_string(t));
      }

      for (int A = 0; A <= maxA; ++A) {
        if (e.YVsA[A] != output.YVsAVsTime()[t][A])
          throw std::runtime_error("Wrong YvsA at " + std::to_string(t));
      }
    }

    for (int i = 0; i < numNuc; ++i) {
      if (nuclib.As()[i] != output.A()[i])
        throw std::runtime_error("Wrong A at nuclide " + std::to_string(i));
      if (nuclib.Zs()[i] != output.Z()[i])
        throw std::runtime_error("Wrong Z at nuclide " + std::to_string(i));
      if (nuclib.BindingEnergiesInMeV()[i] != output.BindingEnergy()[i])
        throw std::runtime_error("Wrong binding energy at nuclide "
            + std::to_string(i));
    }
  }

  // test incremental wirting to disk
  std::size_t bufferSize = NetworkOutput::BufferSize;
  bool exceptionThrow = false;

  auto output = NetworkOutput::CreateNew("out", nuclib, false);

  try {
    // we are going to write 10 * BufferSize + (BufferSize - 1) entries, then
    // abort and check that 10 * BufferSize entries were written to disk

    for (unsigned int t = 0; t < 20 * bufferSize; ++t) {
      if (t >= 10 * bufferSize + bufferSize - 1)
        throw std::runtime_error("simulate an error");

      auto e = MakeEntry(t, numNuc, nuclib.As(), maxA);
      ThermodynamicState thermoState(e.temperature, e.density, e.entropy, 0.0,
          0.0, e.Ye, 0.0);
      output.AddEntry(e.time, { e.dt, "foo" }, 0, 0, thermoState, e.heatingRate,
          e.Y, 1.0);
    }
  } catch (...) {
    exceptionThrow = true;
  }

  if (!exceptionThrow) {
    printf("something went wrong\n");
    return EXIT_FAILURE;
  }

  auto input = NetworkOutput::ReadFromFile("out.h5");

  if (input.NumEntries() != 10 * bufferSize) {
    printf("incremental writing did not work, num entries read: %lu\n",
        input.NumEntries());
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


