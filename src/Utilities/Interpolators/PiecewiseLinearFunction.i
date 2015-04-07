%{
#include "Utilities/Interpolators/InterpolatorBase.hpp"
#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"
%}

%ignore MakeUniquePtr;

%include "Utilities/Interpolators/InterpolatorBase.hpp"
%include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

%template(GeneralPiecewiseLinearFunction_double) GeneralPiecewiseLinearFunction<double>;
%pythoncode %{
  PiecewiseLinearFunction = GeneralPiecewiseLinearFunction_double
%}
