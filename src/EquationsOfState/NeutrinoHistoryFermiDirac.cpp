/// \file NeutrinoHistoryFermiDirac.cpp
/// \author jlippuner
/// \since Feb 4, 2016
///
/// \brief
///
///

#include "EquationsOfState/NeutrinoHistoryFermiDirac.hpp"

void NeutrinoHistoryFermiDirac::PrintInfo(NetworkOutput * const pOutput) const {
  pOutput->Log("# Neutrino History: NeutrinoHistoryThermal\n");
  pOutput->Log("#   Type: %s\n", mIsConst ? "constant" : "time-dependent");
  pOutput->Log("#   Number of distributions: %lu\n", mSpecies.size());

  for (unsigned int i = 0; i < mSpecies.size(); ++i) {
    pOutput->Log("#   Distribution %u:\n", i);
    pOutput->Log("#     Species: %s\n",
        NeutrinoDistribution::NeutrinoSpeciesToStr(mSpecies[i]).c_str());
    pOutput->Log("#     initial temperature: %.5e\n", mT9VsTime.Values()[0][i]);
    pOutput->Log("#     initial eta: %.5e\n", mEtaVsTime.Values()[0][i]);
    pOutput->Log("#     initial normalization: %.5e\n",
        mNormVsTime.Values()[0][i]);
  }
}
