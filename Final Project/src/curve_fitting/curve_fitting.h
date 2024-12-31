#ifndef CURVE_FITTING_H
#define CURVE_FITTING_H

#include <vector>
#include "pp_spline.h"

class CurveFitting {
public:
    // Fit a curve using PPSpline
    static PPSpline fitCurve(const std::vector<double>& x, const std::vector<double>& y, SplineBoundaryCondition bc);

    // Fit a curve on a sphere (spherical curve fitting)
    static PPSpline fitSphericalCurve(const std::vector<std::pair<double, double>>& sphericalPoints, SplineBoundaryCondition bc);
};

#endif // CURVE_FITTING_H
