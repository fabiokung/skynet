%{
#include "Utilities/FunctionVsTime.hpp"
#include "Utilities/Interpolators/InterpolatorBase.hpp"
%}

%include "Utilities/FunctionVsTime.hpp"
%include "Utilities/Interpolators/InterpolatorBase.hpp"

%template(InterpolatorBase_double) InterpolatorBase<double>;
