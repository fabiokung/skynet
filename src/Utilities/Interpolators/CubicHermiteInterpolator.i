%{
#include "Utilities/Interpolators/InterpolatorBase.hpp"
#include "Utilities/Interpolators/CubicHermiteInterpolator.hpp"
%}

%ignore MakeUniquePtr;

%include "Utilities/Interpolators/InterpolatorBase.hpp"
%include "Utilities/Interpolators/CubicHermiteInterpolator.hpp"

%template(CubicHermiteInterpolator_double) CubicHermiteInterpolator<double>;

