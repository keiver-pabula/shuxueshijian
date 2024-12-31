#include "curve_fitting.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fit a curve using PPSpline
PPSpline CurveFitting::fitCurve(const std::vector<double>& x, const std::vector<double>& y, SplineBoundaryCondition bc) {
    if (x.size() != y.size() || x.size() < 2) {
        throw std::invalid_argument("Input data size mismatch or insufficient points.");
    }

    // Convert y to MathFunction for PPSpline
    std::vector<MathFunction> functions;
    for (double yi : y) {
        functions.push_back(MathFunction([yi](double) { return yi; }));
    }

    return PPSpline(1, 3, functions, x.front(), x.back(), x.size(), bc);
}

// Fit a curve on a sphere
PPSpline CurveFitting::fitSphericalCurve(const std::vector<std::pair<double, double>>& sphericalPoints, SplineBoundaryCondition bc) {
    if (sphericalPoints.size() < 2) {
        throw std::invalid_argument("Insufficient spherical points.");
    }

    // Prepare parameter t
    std::vector<double> t(sphericalPoints.size());
    for (size_t i = 0; i < sphericalPoints.size(); ++i) {
        t[i] = static_cast<double>(i);
    }

    // Extract latitude and longitude separately
    std::vector<double> latitudes, longitudes;
    for (const auto& point : sphericalPoints) {
        latitudes.push_back(point.first);
        longitudes.push_back(point.second);
    }

    // Perform cubic spline interpolation for latitudes and longitudes
    PPSpline splineLatitudes = fitCurve(t, latitudes, bc);
    PPSpline splineLongitudes = fitCurve(t, longitudes, bc);

    // Use the coefficients of latitudes and longitudes to construct the spherical spline
    return PPSpline(
        2, // Set dim = 2 for latitude and longitude
        3, // Degree
        {MathFunction([&splineLatitudes](double x) { return splineLatitudes.evaluate(x); }),
         MathFunction([&splineLongitudes](double x) { return splineLongitudes.evaluate(x); })}, // Math functions for lat/lon
        t.front(),
        t.back(),
        t.size(),
        bc
    );
}

