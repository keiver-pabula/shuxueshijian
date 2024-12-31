#include "spline_utils.h"

// Solve a tridiagonal matrix system
std::vector<double> solveTridiagonalMatrix(
    const std::vector<double>& lower,
    const std::vector<double>& diag,
    const std::vector<double>& upper,
    const std::vector<double>& rhs
) {
    int n = diag.size();
    if (n <= 1 || lower.size() != n - 1 || upper.size() != n - 1 || rhs.size() != n) {
        throw std::invalid_argument("Invalid dimensions for tridiagonal matrix.");
    }

    std::vector<double> c_prime(n - 1); // Modified upper diagonal
    std::vector<double> d_prime(n);    // Modified RHS

    // Forward sweep
    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];
    for (int i = 1; i < n; ++i) {
        double m = diag[i] - lower[i - 1] * c_prime[i - 1];
        if (i < n - 1) c_prime[i] = upper[i] / m;
        d_prime[i] = (rhs[i] - lower[i - 1] * d_prime[i - 1]) / m;
    }

    // Back substitution
    std::vector<double> x(n);
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

// Compute the first derivative
double calculateDerivative(double f1, double f2, double x1, double x2) {
    if (x1 == x2) {
        throw std::invalid_argument("Division by zero in derivative calculation.");
    }
    return (f2 - f1) / (x2 - x1);
}

// Compute cumulative chordal length
std::vector<double> computeChordalLength(const std::vector<double>& x, const std::vector<double>& y) {
    int n = x.size();
    if (n != y.size() || n < 2) {
        throw std::invalid_argument("Invalid input size for chordal length calculation.");
    }

    std::vector<double> chordLength(n);
    chordLength[0] = 0.0;

    for (int i = 1; i < n; ++i) {
        double dx = x[i] - x[i - 1];
        double dy = y[i] - y[i - 1];
        chordLength[i] = chordLength[i - 1] + std::sqrt(dx * dx + dy * dy);
    }

    return chordLength;
}

// Generate uniform nodes in an interval
std::vector<double> generateUniformNodes(double start, double end, int numNodes) {
    if (numNodes < 2) {
        throw std::invalid_argument("Number of nodes must be at least 2.");
    }

    std::vector<double> nodes(numNodes);
    double step = (end - start) / (numNodes - 1);

    for (int i = 0; i < numNodes; ++i) {
        nodes[i] = start + i * step;
    }

    return nodes;
}
