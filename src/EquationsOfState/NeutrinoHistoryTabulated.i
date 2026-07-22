%{
#include "Utilities/FunctionVsTime.hpp"
#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "EquationsOfState/NeutrinoDistributionTabulated.hpp"
#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoHistoryTabulated.hpp"
%}

%ignore CreateConstant;
%ignore CreateTimeDependent;
%rename(CreateConstant) CreateConstantSWIG;
%rename(CreateTimeDependent) CreateTimeDependentSWIG;

%include "EquationsOfState/NeutrinoHistoryTabulated.hpp"
