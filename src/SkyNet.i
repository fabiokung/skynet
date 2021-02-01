%module SkyNet

%include <carrays.i>
%include <exception.i>
%include "std_array.i"
%include <std_map.i>
%include <std_string.i>
%include <std_vector.i>
%include <std_shared_ptr.i>
%include "Reactions/ReactionData.i"

// std::array typemaps
%std_array_typemap(double, 24)
%std_array_typemap(double, 7)
%std_array_typemap(int, 2)

%{
#define SWIG_FILE_WITH_INIT

// for exception handler
#include <exception>

// if a class goes into a std::vector, its header must be included here
#include "NuclearData/Nuclide.hpp"
#include "Reactions/REACLIBEntry.hpp"
#include "Reactions/Reaction.hpp"
#include "Reactions/ReactionData.hpp"
#include "Reactions/ReactionLibraryBase.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoHistoryFermiDirac.hpp"
#include "EquationsOfState/NeutrinoHistoryBlackBody.hpp"
#include "Utilities/FunctionVsTime.hpp"

// the SWIG code generates a -Wmaybe-uninitialized warning, so turn it off
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
%}

// instantiate std::vector templates used in SkyNet
%template(std_vector_int) std::vector<int>;
%template(std_vector_double) std::vector<double>;
%template(std_vector_std_vector_double) std::vector<std::vector<double>>;
%template(std_vector_std_string) std::vector<std::string>;
%template(std_vector_const_ReactionLibraryBase_ptr)
    std::vector<const ReactionLibraryBase *>;
%template(std_vector_const_FunctionVsTime_double_ptr)
    std::vector<const FunctionVsTime<double> *>;
%template(std_vector_NeutrinoSpecies) std::vector<NeutrinoSpeciesStruct>;

// we have to ignore the vector constructor from size and resize method for
// classes that don't have default constructors (like Nuclide and Reaction)

%ignore std::vector<Nuclide>::vector(size_type);
%ignore std::vector<Nuclide>::resize(size_type);
%template(std_vector_Nuclide) std::vector<Nuclide>;

%ignore std::vector<REACLIBEntry>::vector(size_type);
%ignore std::vector<REACLIBEntry>::resize(size_type);
%template(std_vector_REACLIBEntry) std::vector<REACLIBEntry>;

%template(ReactionData_Reaction) ReactionData<Reaction>;

%ignore std::vector<Reaction>::vector(size_type);
%ignore std::vector<Reaction>::resize(size_type);
%template(std_vector_Reaction) std::vector<Reaction>;

// instantiate std::array templates used in SkyNet
%template(std_array_double_24) std::array<double, 24>;
%template(std_array_double_7) std::array<double, 7>;
%template(std_array_int_2) std::array<int, 2>;

// instantiate std::map templates used in SkyNet
%template(std_map_std_string_int) std::map<std::string, int>;

%shared_ptr(NeutrinoDistribution)
%shared_ptr(DummyNeutrinoDistribution)
%shared_ptr(NeutrinoDistributionFermiDirac)

%shared_ptr(FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>)
%shared_ptr(NeutrinoHistory)
%shared_ptr(DummyNeutrinoHistory)
%shared_ptr(NeutrinoHistoryFermiDirac)
%shared_ptr(NeutrinoHistoryBlackBody)

%{
// re-enable warnings
#pragma GCC diagnostic pop
%}

// make default exception handler
%exception {
  try {
    $action
  } catch (std::exception& ex) {
    std::string errorMsg = "** SKYNET ERROR ** " + std::string(ex.what());
    SWIG_exception(SWIG_RuntimeError, errorMsg.c_str());
  }
}

// include all interface files
%include "BuildInfo.i"

// FunctionVsTime must appear before it is used in other classes
%include "Utilities/FunctionVsTime.i"

// instantiate FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>
%template(FunctionVsTime_std_shared_ptr_NeutrinoDistribution)
    FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>;

%include "DensityProfiles/ConstantFunction.i"
%include "DensityProfiles/ExpTMinus3.i"
%include "DensityProfiles/PowerLawContinuation.i"
%include "DensityProfiles/Homologous.i"

%include "EquationsOfState/HelmholtzEOS.i"
%include "EquationsOfState/SkyNetScreening.i"
%include "EquationsOfState/NeutrinoDistribution.i"
%include "EquationsOfState/NeutrinoDistributionFermiDirac.i"
%include "EquationsOfState/NeutrinoHistory.i"
%include "EquationsOfState/NeutrinoHistoryFermiDirac.i"
%include "EquationsOfState/NeutrinoHistoryBlackBody.i"

%include "MovieMaker/ColorScale.i"
%include "MovieMaker/ColorScales.i"
%include "MovieMaker/MovieData.i"
%include "MovieMaker/NucleiChart.i"
%include "MovieMaker/NucleiChartOptions.i"
%include "MovieMaker/PiecewiseLinearColorScale.i"
%include "MovieMaker/Scale.i"

%include "Network/NetworkOptions.i"
%include "Network/NetworkOutput.i"
%include "Network/NSE.i"
%include "Network/NSEOptions.i"
%include "Network/ReactionNetwork.i"
%include "Network/TemperatureDensityHistory.i"

%include "NuclearData/Nuclide.i"
%include "NuclearData/NuclideLibrary.i"

%include "Reactions/ConstReactionLibrary.i"
%include "Reactions/REACLIB.i"
%include "Reactions/REACLIBReactionLibrary.i"
%include "Reactions/REACLIBEntry.i"
%include "Reactions/Reaction.i"
%include "Reactions/ReactionPostProcess.i"
%include "Reactions/FFNReactionLibrary.i"
%include "Reactions/Neutrino.i"
%include "Reactions/NeutrinoExternalRatesReactionLibrary.i"
%include "Reactions/NeutrinoReactionLibrary.i"

%include "Utilities/Constants.i"
%include "Utilities/Interpolators/InterpolatorBase.i"
%include "Utilities/Interpolators/PiecewiseLinearFunction.i"
%include "Utilities/Interpolators/CubicHermiteInterpolator.i"
