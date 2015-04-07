%{
#include "Utilities/FunctionVsTime.hpp"
%}

%ignore MakeUniquePtr;

%include "Utilities/FunctionVsTime.hpp"
%template(FunctionVsTime_double) FunctionVsTime<double>;
