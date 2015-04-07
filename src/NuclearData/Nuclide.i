#pragma SWIG nowarn=305
%warnfilter(302) Nuclide::PartitionFunctionT9_SWIG;
%warnfilter(302) Nuclide::PartitionFunctionT9;

%{
#include "NuclearData/Nuclide.hpp"
%}

%ignore NumPartitionFunctionEntries;

%include "NuclearData/Nuclide.hpp"

