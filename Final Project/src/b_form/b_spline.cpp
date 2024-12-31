#include "b_spline.h"
#include <stdexcept>
#include <cmath>
#include <vector>
#include <iostream>
#include <ostream>
#include <iomanip>
#include "../pp_form/piecewise_polynomial.h"

// BSpline::BSpline(const std::vector<double>& controlPoints, int degree)
//     : controlPoints(controlPoints), degree(degree) {
//     if (controlPoints.size() < degree + 1) {
//         throw std::invalid_argument("Not enough control points for the given degree.");
//     }
//     computeKnots(controlPoints); 
// }

// Constructor for Task B
BSpline::BSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_points, SplineBoundaryCondition bc)
    : degree(order - 1) {
    // Generate evenly spaced points
    std::vector<double> x_values(num_points);
    for (int i = 0; i < num_points; ++i) {
        x_values[i] = a + i * (b - a) / (num_points - 1);
    }

    // Compute control points using the provided functions
    controlPoints.resize(num_points);
    for (int i = 0; i < num_points; ++i) {
        controlPoints[i] = functions[0].evaluate(x_values[i]);
    }

    // Determine the knot vector based on the boundary condition
    switch (bc) {
        case NATURAL_SPLINE:
            computeNaturalKnots(x_values);
            break;
        case CLAMPED:
            computeClampedKnots(x_values);
            break;
        case PERIODIC_CONDITION:
            computePeriodicKnots(x_values);
            break;
        default:
            throw std::invalid_argument("Unsupported boundary condition");
    }
}

BSpline::BSpline(const std::vector<double>& controlPoints, int degree)
    : controlPoints(controlPoints), degree(degree) {
    if (controlPoints.size() < degree + 1) {
        throw std::invalid_argument("Not enough control points for the given degree.");
    }
    computeKnots(controlPoints);
}

// Constructor with explicit knots
BSpline::BSpline(const std::vector<double>& controlPoints, const std::vector<double>& knots, int degree)
    : controlPoints(controlPoints), knots(knots), degree(degree) {
    if (controlPoints.size() + degree + 1 != knots.size()) {
        throw std::invalid_argument("Invalid number of knots for the given degree and control points.");
    }
}


void BSpline::computePeriodicKnots(const std::vector<double>& x) {
    int numKnots = x.size() + degree + 1;
    knots.resize(numKnots);

    for (int i = 0; i <= degree; ++i) {
        knots[i] = x[x.size() - degree + i];  // Periodic start knots
        knots[numKnots - 1 - i] = x[i];  // Periodic end knots
    }

    for (int i = 0; i < x.size(); ++i) {
        knots[degree + i] = x[i];
    }

    std::cout << "Periodic Knot Vector: ";
    for (double knot : knots) {
        std::cout << knot << " ";
    }
    std::cout << "\n";
}

void BSpline::computeClampedKnots(const std::vector<double>& x) {
    int numKnots = x.size() + degree + 1;
    knots.resize(numKnots);

    for (int i = 0; i <= degree; ++i) {
        knots[i] = x.front();  // Clamped start knots
        knots[numKnots - 1 - i] = x.back();  // Clamped end knots
    }

    for (int i = 1; i < x.size() - degree; ++i) {
        knots[degree + i] = x[i];
    }

    std::cout << "Clamped Knot Vector: ";
    for (double knot : knots) {
        std::cout << knot << " ";
    }
    std::cout << "\n";
}

void BSpline::computeNaturalKnots(const std::vector<double>& x) {
    int numKnots = x.size() + degree + 1;
    knots.resize(numKnots);

    for (int i = 0; i <= degree; ++i) {
        knots[i] = x.front();  // Start knots
        knots[numKnots - 1 - i] = x.back();  // End knots
    }

    double interval = (x.back() - x.front()) / (x.size() - 1);
    for (int i = 1; i < x.size() - degree; ++i) {
        knots[degree + i] = knots[degree] + i * interval;
    }

    std::cout << "Natural Knot Vector: ";
    for (double knot : knots) {
        std::cout << knot << " ";
    }
    std::cout << "\n";
}

std::vector<double> BSpline::fitCurve(const std::vector<double>& x, const std::string& method) {
    std::vector<double> parameters;

    if (method == "uniform") {
        double step = (x.back() - x.front()) / (x.size() - 1);
        for (size_t i = 0; i < x.size(); ++i) {
            parameters.push_back(x.front() + i * step);
        }
    } else if (method == "chordal") {
        parameters.push_back(0.0);
        for (size_t i = 1; i < x.size(); ++i) {
            parameters.push_back(parameters.back() + std::abs(x[i] - x[i - 1]));
        }
        for (double& p : parameters) {
            p /= parameters.back();
        }
    } else {
        throw std::invalid_argument("Unknown parameterization method: " + method);
    }

    return parameters;
}

void BSpline::computeKnots(const std::vector<double>& x) {
    int numKnots = x.size() + degree + 1;
    knots.resize(numKnots);

    // For Curve 1, use chordal parameterization
    if (!x.empty() && x.front() == 0.0 && x.back() == 1.0) { 
        for (int i = 0; i <= degree; ++i) {
            knots[i] = x.front();
            knots[numKnots - 1 - i] = x.back();
        }
        for (int i = 1; i < x.size() - degree; ++i) {
            knots[degree + i] = x[i];
        }
    } else {
        // Default uniform knot vector
        for (int i = 0; i <= degree; ++i) {
            knots[i] = x.front();
            knots[numKnots - 1 - i] = x.back();
        }
        for (int i = 1; i < x.size() - degree; ++i) {
            knots[degree + i] = x.front() + i * (x.back() - x.front()) / (x.size() - 1);
        }
    }

    std::cout << "Updated Knot Vector: ";
    for (double knot : knots) {
        std::cout << knot << " ";
    }
    std::cout << "\n";
}


bool debugBasis = false; // Set to true if debugging is needed

double BSpline::evaluateBasis(int i, int k, double t) const {
    if (k == 0) {
        double value = (t >= knots[i] && t < knots[i + 1]) ? 1.0 : 0.0;
        if (debugBasis) {
            std::cout << "Basis (i=" << i << ", k=" << k << ", t=" << t << "): " << value << std::endl;
        }
        return value;
    }

    double left = 0.0, right = 0.0;
    if (knots[i + k] != knots[i]) {
        left = (t - knots[i]) / (knots[i + k] - knots[i]) * evaluateBasis(i, k - 1, t);
    }
    if (knots[i + k + 1] != knots[i + 1]) {
        right = (knots[i + k + 1] - t) / (knots[i + k + 1] - knots[i + 1]) * evaluateBasis(i + 1, k - 1, t);
    }

    double value = left + right;
    if (debugBasis) {
        std::cout << "Basis (i=" << i << ", k=" << k << ", t=" << t << "): " << value << std::endl;
    }
    return value;
}

double BSpline::evaluate(double t) const {
    if (t < knots[degree] || t > knots[knots.size() - degree - 1]) {
        throw std::out_of_range("t is out of range of the B-spline.");
    }

    double result = 0.0;
    std::cout << "Evaluating B-Spline at t = " << t << std::endl;
    for (size_t i = 0; i < controlPoints.size(); ++i) {
        double basis = evaluateBasis(i, degree, t);
        result += controlPoints[i] * basis;        
    }
    std::cout << "Result: " << result << std::endl;
    return result;
}

void BSpline::fit(const std::vector<double>& x, const std::vector<double>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y must have the same size.");
    }

    controlPoints = y; // Temporarily set control points to y-values
    computeKnots(x);
}


void BSpline::print(std::ostream& os) const {
    os << std::fixed << std::setprecision(6); // Consistent decimal precision

    // Print knots
    os << "Knots: ";
    for (const auto& knot : knots) {
        os << knot << " ";
    }
    os << "\n";

    // Print piecewise polynomials
    os << "Piecewise Polynomials:\n";

    // Loop over intervals
    for (size_t i = degree; i < knots.size() - 1 - degree; ++i) {
        // Define the interval
        double interval_start = knots[i];
        double interval_end = knots[i + 1];
        os << "Interval: [" << interval_start << ", " << interval_end << "]\n";

        // Compute the polynomial coefficients for this interval
        std::vector<double> coeffs(degree + 1, 0.0);
        for (int j = 0; j <= degree; ++j) {
            // Accumulate contributions from each control point
            coeffs[j] = controlPoints[j] * evaluateBasis(i - degree + j, degree, (interval_start + interval_end) / 2.0);
        }

        os << "Polynomial Coefficients: ";
        for (size_t j = 0; j < coeffs.size(); ++j) {
            os << coeffs[j] << (j < coeffs.size() - 1 ? ", " : "\n");
        }
    }
}

void BSpline::printDetailed(std::ostream& os) const {
    os << std::fixed << std::setprecision(6); // Consistent decimal precision

    // Print knots
    os << "Knots: ";
    for (const auto& knot : knots) {
        os << knot << " ";
    }
    os << "\n";

    // Print piecewise polynomials
    os << "Piecewise Polynomials:\n";

    // Loop over intervals
    for (size_t i = degree; i < knots.size() - 1 - degree; ++i) {
        // Define the interval
        double interval_start = knots[i];
        double interval_end = knots[i + 1];
        os << "Interval: [" << interval_start << ", " << interval_end << "]\n";

        // Compute the polynomial coefficients for this interval
        std::vector<double> coeffs(degree + 1, 0.0);
        for (int j = 0; j <= degree; ++j) {
            coeffs[j] = 0.0; // Clear previous coefficient calculation
            for (size_t k = 0; k < controlPoints.size(); ++k) {
                double basis = evaluateBasis(k, j, (interval_start + interval_end) / 2.0);
                coeffs[j] += controlPoints[k] * basis;
            }
        }

        os << "Polynomial Coefficients: ";
        for (size_t j = 0; j < coeffs.size(); ++j) {
            os << coeffs[j] << (j < coeffs.size() - 1 ? ", " : "\n");
        }
    }
}

void BSpline::printDebug(std::ostream& os) const {
    // Print knots
    os << "Knots:\n";
    for (const auto& knot : knots) {
        os << knot << " ";
    }
    os << "\n";

    os << "Control Points:\n";
    for (const auto& cp : controlPoints) {
        os << cp << " ";
    }
    os << "\n";

    // Check polynomials if the `segments` exist
    os << "\nPiecewise Polynomials:\n";
    if (segments.empty()) {
        os << "  [Error] Polynomial segments are not computed.\n";
    } else {
        for (size_t i = 0; i < segments.size(); ++i) {
            os << "  Segment " << i + 1 << ":\n";
            const auto& coeffs = segments[i].getValues(); // Example accessor; replace with your implementation
            os << "    Coefficients: ";
            for (const auto& coeff : coeffs) {
                os << coeff << " ";
            }
            os << "\n";
        }
    }
}
