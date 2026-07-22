add_subdirectory(Helmholtz)

set(EquationsOfState_SOURCES
  HelmholtzEOS.cpp
  NeutrinoDistributionFermiDirac.cpp
  NeutrinoHistoryFermiDirac.cpp
  NeutrinoHistoryBlackBody.cpp
  TabulatedNeutrinoSpectrum.cpp
  NeutrinoDistributionTabulated.cpp
  NeutrinoHistoryTabulated.cpp
  SkyNetScreening.cpp
)

add_SkyNet_library(EquationsOfState "${EquationsOfState_SOURCES}")

set(EquationsOfState_SWIG_DEPS
  EOS.hpp
  HelmholtzEOS.i
  HelmholtzEOS.hpp
  Screening.hpp
  SkyNetScreening.hpp
  SkyNetScreening.i
  NeutrinoDistribution.i
  NeutrinoDistribution.hpp
  NeutrinoHistory.i
  NeutrinoHistory.hpp
  NeutrinoDistributionFermiDirac.i
  NeutrinoDistributionFermiDirac.hpp
  NeutrinoHistoryFermiDirac.i
  NeutrinoHistoryFermiDirac.hpp
  NeutrinoHistoryBlackBody.i
  NeutrinoHistoryBlackBody.hpp
  TabulatedNeutrinoSpectrum.hpp
  NeutrinoDistributionTabulated.i
  NeutrinoDistributionTabulated.hpp
  NeutrinoHistoryTabulated.i
  NeutrinoHistoryTabulated.hpp
)

add_SWIG_dependencies("${EquationsOfState_SWIG_DEPS}")
