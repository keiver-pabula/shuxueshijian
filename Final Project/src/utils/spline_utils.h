#ifndef SPLINE_UTILS_H
#define SPLINE_UTILS_H

#include <vector>
#include <stdexcept>
#include <cmath>

// Solve a tridiagonal matrix system using the Thomas algorithm
std::vector<double> solveTridiagonalMatrix(
    const std::vector<double>& lower,  // Sub-diagonal elements
    const std::vector<double>& diag,   // Main diagonal elements
    const std::vector<double>& upper,  // Super-diagonal elements
    const std::vector<double>& rhs     // Right-hand side
);

// Compute the first derivative
double calculateDerivative(double f1, double f2, double x1, double x2);

// Compute cumulative chordal length
std::vector<double> computeChordalLength(const std::vector<double>& x, const std::vector<double>& y);

// Generate uniform nodes in an interval
std::vector<double> generateUniformNodes(double start, double end, int numNodes);

#endif // SPLINE_UTILS_H
