/// \file Nuclide.hpp
/// \author jlippuner
/// \since Jul 2, 2014
///
/// \brief
///
///

#ifndef SKYNET_NUCLEARDATA_NUCLIDE_HPP_
#define SKYNET_NUCLEARDATA_NUCLIDE_HPP_

#include <array>
#include <string>
#include <vector>

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

#define NUCLIDE_PARTITION_FUNCTION_T9                                           \
0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,                              \
1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0

class Nuclide {
public:
  Nuclide(const int z, const int a, const double massExcess,
      const double groundStateSpin,
      const std::array<double, 24>& partitionFunctionLog10,
      const std::string& name);

  static std::string GetName(const int Z, const int A);

  static std::array<int, 2> GetAAndZ(const std::string& name);

  static Nuclide CreateDummy(const int A, const std::string& name);

  // SWIG does not know how to convert double[24] to a Python sequence, so
  // we provide a std::array for SWIG, this must appear before
  // PartitionFunctionT9
#ifdef SWIG
%rename(PartitionFunctionT9) PartitionFunctionT9_SWIG;
#endif // SWIG
  static constexpr std::array<double, 24> PartitionFunctionT9_SWIG {{
    NUCLIDE_PARTITION_FUNCTION_T9
  }};

  // Intel cannot handle static constexpr std::array, so we have to do this the
  // old fashioned way
  static constexpr double PartitionFunctionT9[24] = {
      NUCLIDE_PARTITION_FUNCTION_T9
  };

  static constexpr std::size_t NumPartitionFunctionEntries = 24;

  // The maximum number of the log10 of the partition function, where it is
  // capped out in the input data
  static constexpr double PartitionFunctionLog10Cap = 5.0;

  int Z() const {
    return mZ;
  }

  int A() const {
    return mA;
  }

  const std::string& Name() const {
    return mName;
  }

  double MassExcess() const {
    return mMassExcess;
  }

  double GroundStateSpin() const {
    return mGroundStateSpin;
  }

  const std::array<double, 24>& PartitionFunctionLog10() const {
    return mPartitionFunctionLog10;
  }

private:
  int mZ;
  int mA;
  std::string mName;
  double mMassExcess; // in MeV
  double mGroundStateSpin;
  std::array<double, 24> mPartitionFunctionLog10;
};

#endif // SKYNET_NUCLEARDATA_NUCLIDE_HPP_
