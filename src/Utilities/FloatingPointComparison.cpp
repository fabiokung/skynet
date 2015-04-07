
/// \file FloatingPointComparison.cpp
/// \author jlippuner
/// \since Jan 16, 2014
///
/// \brief Implement floating point comparison functions.
///
///

#include "FloatingPointComparison.hpp"

// copied from gsl/sys/fcmp.c so that the gsl library is not required to compile
int skynet_gsl_fcmp(const double x1, const double x2, const double epsilon) {
    int exponent;
    double delta, difference;

    /* Find exponent of largest absolute value */
    {
        double max = (fabs(x1) > fabs(x2)) ? x1 : x2;
        frexp(max, &exponent);
    }

    /* Form a neighborhood of size  2 * delta */
    delta = ldexp(epsilon, exponent);
    difference = x1 - x2;

    if (difference > delta) /* x1 > x2 */
        return 1;
    else if (difference < -delta) /* x1 < x2 */
        return -1;
    else /* -delta <= difference <= delta */
        return 0; /* x1 ~=~ x2 */
}

/// function to compare doubles, returns:
///  0 if x = y (within tolerance)
/// +1 if x > y
/// -1 if x < y
int fcmp(const double x, const double y, const double epsilon) {
    // need to treat 0.0 differently, because gsl_fcmp cannot handle it
    if (y == 0.0) {
        if (x == 0.0)
            return 0;
        else
            return -fcmp(y, x, epsilon);
    }
    if (x == 0.0) {
        if (y > epsilon)
            return -1;
        if (y < -epsilon)
            return 1;
        return 0;
    }
    return skynet_gsl_fcmp(x, y, epsilon);
}
