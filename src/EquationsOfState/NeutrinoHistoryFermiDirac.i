%{
#include "Utilities/FunctionVsTime.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoHistoryFermiDirac.hpp"
%}

%include "Utilities/FunctionVsTime.hpp"
%include "EquationsOfState/NeutrinoDistribution.hpp"
%include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"
%include "EquationsOfState/NeutrinoHistory.hpp"

%ignore CreateConstant;
%ignore CreateTimeDependent;
%rename(CreateConstant) CreateConstantSWIG;
%rename(CreateTimeDependent) CreateTimeDependentSWIG;

%include "EquationsOfState/NeutrinoHistoryFermiDirac.hpp"
