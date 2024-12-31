#include "cubic_spline.h"
#include "../utils/spline_utils.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <ostream>

void CubicSpline::solveTridiagonal(const std::vector<double>& h, const std::vector<double>& alpha) {
    int n = x.size() - 1;
    std::vector<double> lower(n), diag(n + 1), upper(n);
    std::vector<double> rhs(n + 1);

    // Fill the tridiagonal matrix
    for (int i = 1; i < n; ++i) {
        lower[i - 1] = h[i - 1];
        diag[i] = 2 * (h[i - 1] + h[i]);
        upper[i - 1] = h[i];
        rhs[i] = alpha[i];
    }

    // Natural boundary conditions
    diag[0] = 1.0;
    diag[n] = 1.0;
    rhs[0] = 0.0;
    rhs[n] = 0.0;

    // Solve for c
    c = solveTridiagonalMatrix(lower, diag, upper, rhs);

    // Solve for b and d
    b.resize(n);
    d.resize(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
}

void CubicSpline::fitNatural(const std::vector<double>& x_points, const std::vector<double>& y_points) {
    x = x_points;
    y = y_points;
    int n = x.size() - 1;

    std::vector<double> h(n), alpha(n - 1);
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    for (int i = 1; i < n; ++i) {
        alpha[i - 1] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }

    // Natural boundary conditions
    std::vector<double> lower(n), diag(n + 1), upper(n), rhs(n + 1);
    diag[0] = 1.0;          // Boundary condition at the start
    diag[n] = 1.0;          // Boundary condition at the end
    rhs[0] = 0.0;           // Second derivative at the start
    rhs[n] = 0.0;           // Second derivative at the end

    for (int i = 1; i < n; ++i) {
        lower[i - 1] = h[i - 1];
        diag[i] = 2 * (h[i - 1] + h[i]);
        upper[i - 1] = h[i];
        rhs[i] = alpha[i - 1];
    }

    c = solveTridiagonalMatrix(lower, diag, upper, rhs);

    // Solve for b and d
    b.resize(n);
    d.resize(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }

    std::cout << "Cubic Spline Coefficients:\n";
    for (size_t i = 0; i < a.size(); ++i) {
        std::cout << "a[" << i << "] = " << a[i]
                  << ", b[" << i << "] = " << b[i]
                  << ", c[" << i << "] = " << c[i]
                  << ", d[" << i << "] = " << d[i] << "\n";
    }
}

void CubicSpline::fitClamped(const std::vector<double>& x_points, const std::vector<double>& y_points, double df0, double dfn) {
    x = x_points;
    y = y_points;
    int n = x.size() - 1;

    std::vector<double> h(n), alpha(n + 1);
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    alpha[0] = 3 * ((y[1] - y[0]) / h[0] - df0);
    alpha[n] = 3 * (dfn - (y[n] - y[n - 1]) / h[n - 1]);
    for (int i = 1; i < n; ++i) {
        alpha[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }

    // Clamped boundary conditions
    std::vector<double> lower(n), diag(n + 1), upper(n), rhs(n + 1);
    diag[0] = 2 * h[0];
    diag[n] = 2 * h[n - 1];
    rhs[0] = alpha[0];
    rhs[n] = alpha[n];

    for (int i = 1; i < n; ++i) {
        lower[i - 1] = h[i - 1];
        diag[i] = 2 * (h[i - 1] + h[i]);
        upper[i - 1] = h[i];
        rhs[i] = alpha[i];
    }

    c = solveTridiagonalMatrix(lower, diag, upper, rhs);

    b.resize(n);
    d.resize(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
}

void CubicSpline::fitPeriodic(const std::vector<double>& x_points, const std::vector<double>& y_points) {
    x = x_points;
    y = y_points;
    int n = x.size() - 1;

    std::vector<double> h(n), alpha(n + 1);
    for (int i = 0; i < n; ++i) {
        h[i] = x[i + 1] - x[i];
    }
    for (int i = 1; i < n; ++i) {
        alpha[i] = 3 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
    }
    alpha[0] = 3 * ((y[1] - y[0]) / h[0] - (y[n] - y[n - 1]) / h[n - 1]);

    // Periodic boundary conditions
    std::vector<double> lower(n), diag(n + 1), upper(n), rhs(n + 1);
    diag[0] = 2 * h[0];
    diag[n] = 2 * h[n - 1];
    rhs[0] = alpha[0];
    rhs[n] = alpha[n];

    for (int i = 1; i < n; ++i) {
        lower[i - 1] = h[i - 1];
        diag[i] = 2 * (h[i - 1] + h[i]);
        upper[i - 1] = h[i];
        rhs[i] = alpha[i];
    }

    c = solveTridiagonalMatrix(lower, diag, upper, rhs);

    b.resize(n);
    d.resize(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
    }
}

double CubicSpline::evaluate(double t) const {
    if (t < x.front()) {
        // Return the value at the first segment
        double dx = t - x.front();
        return y.front() + b[0] * dx + c[0] * dx * dx + d[0] * dx * dx * dx;
    } else if (t > x.back()) {
        // Return the value at the last segment
        double dx = t - x.back();
        size_t n = x.size() - 1;
        return y[n] + b[n - 1] * dx + c[n - 1] * dx * dx + d[n - 1] * dx * dx * dx;
    }

    // Search for the correct segment
    for (size_t i = 0; i < x.size() - 1; ++i) {
        if (t >= x[i] && t <= x[i + 1]) {
            double dx = t - x[i];
            return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    }

    throw std::out_of_range("Value out of bounds.");
}



void CubicSpline::fitCurve(const std::vector<double>& x_points, const std::vector<double>& y_points, bool useChordalLength) {
    if (x_points.size() != y_points.size() || x_points.size() < 2) {
        throw std::invalid_argument("Invalid input sizes for curve fitting.");
    }

    // Compute parameters based on chordal length or uniform distribution
    std::vector<double> parameters;
    if (useChordalLength) {
        parameters = computeChordalLength(x_points, y_points); // Use chordal length
    } else {
        // Scale uniform nodes to match the range of x_points
        double start = x_points.front();
        double end = x_points.back();
        parameters = generateUniformNodes(start, end, x_points.size()); // Map uniform nodes to x range
    }

    // Fit the spline using the calculated parameters
    fitNatural(parameters, y_points); // Use natural cubic spline for fitting
}
