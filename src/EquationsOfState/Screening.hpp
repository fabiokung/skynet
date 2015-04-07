/// \file Screening.hpp
/// \author jlippuner
/// \since Feb 22, 2017
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_SCREENING_HPP_
#define SRC_EQUATIONSOFSTATE_SCREENING_HPP_

#include <vector>

#include "EquationsOfState/EOS.hpp"
#include "Network/NetworkOutput.hpp"
#include "NuclearData/NuclideLibrary.hpp"

class ScreeningState {
public:
  ScreeningState(const std::vector<double>& mus) :
    mMus(mus) {}

  virtual ~ScreeningState() {}

  inline double MuScreening(const int Z) const {
    return mMus[Z];
  }

  const std::vector<double>& Mus() const {
    return mMus;
  }

private:
  std::vector<double> mMus;
};

class Screening {
public:
  Screening(const NuclideLibrary& nuclib) :
    mMaxZ(0) {
    for (int i = 0; i < nuclib.NumNuclides(); ++i) {
      mMaxZ = std::max(mMaxZ, nuclib.Zs()[i]);
    }
  }

  Screening(const int maxZ) :
    mMaxZ(maxZ) {}

  virtual ~Screening() {}

  virtual std::unique_ptr<Screening> MakeUniquePtr() const =0;

  virtual void PrintInfo(NetworkOutput * const pOutput) const =0;

  int MaxZ() const {
    return mMaxZ;
  }

  void SetMaxZ(const int maxZ) {
    mMaxZ = maxZ;
  }

  virtual ScreeningState Compute(const ThermodynamicState thermoState,
      const std::vector<double>& Y,
      const NuclideLibrary& nuclideLibrary) const =0;

protected:
  int mMaxZ; // the maximum Z for which screening corrections will be computed
};

#endif // SRC_EQUATIONSOFSTATE_SCREENING_HPP_
