/// \file SkyNetScreening.hpp
/// \author jlippuner
/// \since Feb 22, 2017
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_SKYNETSCREENING_HPP_
#define SRC_EQUATIONSOFSTATE_SKYNETSCREENING_HPP_

#include "EquationsOfState/Screening.hpp"

class SkyNetScreening : public Screening {
public:
  SkyNetScreening(const NuclideLibrary& nuclib) :
    Screening(nuclib) {}

  SkyNetScreening(const int maxZ) :
    Screening(maxZ) {}

  virtual std::unique_ptr<Screening> MakeUniquePtr() const {
    return std::unique_ptr<Screening>(new SkyNetScreening(*this));
  }

  virtual void PrintInfo(NetworkOutput * const pOutput) const;

  virtual ScreeningState Compute(const ThermodynamicState thermoState,
      const std::vector<double>& Y, const NuclideLibrary& nuclideLibrary) const;
};

#endif // SRC_EQUATIONSOFSTATE_SKYNETSCREENING_HPP_
