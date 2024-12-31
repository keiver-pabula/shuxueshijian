#include "pp_spline.h"
#include <stdexcept>
#include <iomanip>
#include <vector>

PPSpline::PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree, SplineBoundaryCondition bc)
    : intervals(intervals), degree(degree) {
    if (intervals.size() != values.size() + 1) {
        throw std::invalid_argument("Invalid number of intervals and values.");
    }
    if (degree != 1 && degree != 3) {
        throw std::invalid_argument("Only linear (degree=1) and cubic (degree=3) splines are supported.");
    }

    if (degree == 1) {
        computeLinearSpline(values,0);
    } else if (degree == 3) {
        computeCubicSpline(values,0, bc);
    }
}

PPSpline::PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree, SplineBoundaryCondition bc, bool equalSizes)
    : intervals(intervals), degree(degree) {
    if (!equalSizes && intervals.size() != values.size() + 1) {
        throw std::invalid_argument("Invalid number of intervals and values.");
    }

    if (degree != 1 && degree != 3) {
        throw std::invalid_argument("Only linear (degree=1) and cubic (degree=3) splines are supported.");
    }

    // If sizes are equal, skip interval adjustment
    std::vector<double> adjustedIntervals = intervals;
    if (!equalSizes) {
        adjustedIntervals.pop_back();
    }

    if (degree == 1) {
        computeLinearSpline(values, 0);
    } else if (degree == 3) {
        computeCubicSpline(values, 0, bc);
    }
}


PPSpline::PPSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_intervals, SplineBoundaryCondition bc)
    : degree(order) {
    std::cout << "PPSpline constructor: dim = " << dim << ", functions.size() = " << functions.size() << std::endl;
    if (functions.size() != static_cast<size_t>(dim)) {
        throw std::invalid_argument("Number of functions must match the dimension.");
    }

    // Generate intervals
    double step = (b - a) / (num_intervals - 1);
    for (int i = 0; i <= num_intervals; ++i) {
        intervals.push_back(a + i * step);
    }

    // Compute values
    std::vector<double> values(num_intervals + 1);
    for (size_t i = 0; i < intervals.size(); ++i) {
        values[i] = functions[0].evaluate(intervals[i]);
    }

    // Compute spline based on the order
    if (order == 1) {
        computeLinearSpline(values,0);
    } else if (order == 3) {
        computeCubicSpline(values,0, bc);
    } else {
        throw std::invalid_argument("Unsupported spline order.");
    }
}

PPSpline::PPSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_intervals, SplineBoundaryCondition bc, bool is_point_2)
    : degree(order) {
    std::cout << "PPSpline constructor: dim = " << dim << ", functions.size() = " << functions.size() << std::endl;
    if (functions.size() != static_cast<size_t>(dim)) {
        throw std::invalid_argument("Number of functions must match the dimension.");
    }

    // Generate intervals
    double step = (b - a) / (num_intervals - 1);
    for (int i = 0; i <= num_intervals; ++i) {
        intervals.push_back(a + i * step);
    }

    // Resize coefficients to match the number of dimensions
    coefficients.resize(dim);

    // Compute splines for each dimension
    for (int d = 0; d < dim; ++d) {
        std::vector<double> values(num_intervals + 1);
        for (size_t i = 0; i < intervals.size(); ++i) {
            values[i] = functions[d].evaluate(intervals[i]);
        }

        if (order == 1) {
            computeLinearSpline(values, d);
        } else if (order == 3) {
            computeCubicSpline(values, d, bc);
        } else {
            throw std::invalid_argument("Unsupported spline order.");
        }
    }
}


PPSpline::PPSpline(int dim, int order, const std::vector<MathFunction>& functions, 
                   const std::vector<double>& time_points, SplineBoundaryCondition bc, 
                   double da, double db)
    : degree(order) {
    if (functions.size() != static_cast<size_t>(dim)) {
        throw std::invalid_argument("Number of functions must match the dimension.");
    }
    if (time_points.size() < 2) {
        throw std::invalid_argument("Invalid number of time points.");
    }

    intervals = time_points; // Use provided time points as intervals

    // Compute values of the function at each time point
    std::vector<double> values(time_points.size());
    for (size_t i = 0; i < time_points.size(); ++i) {
        values[i] = functions[0].evaluate(time_points[i]);
    }

    // Compute the cubic spline based on the boundary condition
    if (degree == 3) {
        computeCubicSplineForPoint2(values, bc,0.0,0.0);
    } else {
        throw std::invalid_argument("Only cubic splines are supported for this constructor.");
    }
}

PPSpline::PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree, SplineBoundaryCondition bc, double da, double db)
    : intervals(intervals), degree(degree), da(da), db(db) {
    if (intervals.size() != values.size() + 1) {
        throw std::invalid_argument("Invalid number of intervals and values.");
    }
    if (degree != 1 && degree != 3) {
        throw std::invalid_argument("Only linear (degree=1) and cubic (degree=3) splines are supported.");
    }

    if (degree == 1) {
        computeLinearSpline(values, 0);
    } else if (degree == 3) {
        computeCubicSpline(values,0 , bc);
    }
}

void PPSpline::computeLinearSpline(const std::vector<double>& values, int dim) {
    for (size_t i = 0; i < values.size() - 1; ++i) {
        double slope = (values[i + 1] - values[i]) / (intervals[i + 1] - intervals[i]);
        coefficients.push_back(slope);
        coefficients.push_back(values[i] - slope * intervals[i]);
    }
}

void PPSpline::computeCubicSpline(const std::vector<double>& values, int dim, SplineBoundaryCondition bc) {
    size_t n = values.size();
    std::vector<double> h(n - 1);
    std::vector<double> alpha(n - 1);

    // Step 1: Calculate intervals (h) and alpha
    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = intervals[i + 1] - intervals[i];
        if (i > 0) {
            alpha[i] = 3 * (values[i + 1] - values[i]) / h[i] - 3 * (values[i] - values[i - 1]) / h[i - 1];
        }
    }

    // Step 2: Solve tridiagonal system
    std::vector<double> l(n, 1.0), mu(n, 0.0), z(n, 0.0);
    if (bc == NATURAL_SPLINE) {
        l[0] = 1.0;
        l[n - 1] = 1.0;
    } else if (bc == CLAMPED) {
        l[0] = 2 * h[0];
        mu[0] = 0.5;
        z[0] = 3 * (values[1] - values[0]) / h[0];
        l[n - 1] = 2 * h[n - 2];
        z[n - 1] = 3 * (values[n - 1] - values[n - 2]) / h[n - 2];
    }

    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2 * (intervals[i + 1] - intervals[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    // Step 3: Compute coefficients
    std::vector<double> c(n, 0.0), b(n - 1), d(n - 1);
    for (int j = n - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (values[j + 1] - values[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    // Store coefficients
    for (size_t i = 0; i < n - 1; ++i) {
        coefficients.push_back(values[i]);  // a
        coefficients.push_back(b[i]);      // b
        coefficients.push_back(c[i]);      // c
        coefficients.push_back(d[i]);      // d
    }
}

double PPSpline::evaluate(double x) const {
    for (size_t i = 0; i < intervals.size() - 1; ++i) {
        if (x >= intervals[i] && x < intervals[i + 1]) {
            if (degree == 1) {
                // Linear evaluation
                double slope = coefficients[2 * i];
                double intercept = coefficients[2 * i + 1];
                return slope * x + intercept;
            } else if (degree == 3) {
                // Cubic evaluation
                double dx = x - intervals[i];
                double a = coefficients[4 * i];
                double b = coefficients[4 * i + 1];
                double c = coefficients[4 * i + 2];
                double d = coefficients[4 * i + 3];
                return a + b * dx + c * dx * dx + d * dx * dx * dx;
            }
        }
    }
    throw std::out_of_range("x is outside the spline interval.");
}

void PPSpline::print(std::ostream& os) const {
    os << std::fixed << std::setprecision(6);
    os << "Piecewise-Polynomial Spline (pp-Form):\n";
    for (size_t i = 0; i < intervals.size() - 1; ++i) {
        os << "Interval: [" << intervals[i] << ", " << intervals[i + 1] << "]\n";
        if (degree == 1) {
            os << "y = " << coefficients[2 * i] << "x + " << coefficients[2 * i + 1] << "\n";
        } else if (degree == 3) {
            os << "y = " << coefficients[4 * i] << " + " << coefficients[4 * i + 1] << "x + "
               << coefficients[4 * i + 2] << "x^2 + " << coefficients[4 * i + 3] << "x^3\n";
        }
    }
}

void PPSpline::printDetailed(std::ostream& os) const {
    piecewisePolynomial.print(os);
}


void PPSpline::computeCubicSplineForPoint2(const std::vector<double>& values, SplineBoundaryCondition bc, double da, double db) {
    size_t n = values.size();
    std::vector<double> h(n - 1), alpha(n - 1);

    // Step 1: Calculate first and second divided differences
    for (size_t i = 0; i < n - 1; ++i) {
        h[i] = intervals[i + 1] - intervals[i];
        if (i > 0) {
            alpha[i] = 3 * (values[i + 1] - values[i]) / h[i] - 3 * (values[i] - values[i - 1]) / h[i - 1];
        }
    }

    // Step 2: Set up tridiagonal system
    std::vector<double> l(n, 1.0), mu(n - 1), z(n, 0.0);

    if (bc == NATURAL_SPLINE) {
        l[0] = 1.0;
        l[n - 1] = 1.0;
    } else if (bc == CLAMPED) {
        l[0] = 2 * h[0];
        z[0] = 3 * (values[1] - values[0]) / h[0] - da;
        l[n - 1] = 2 * h[n - 2];
        z[n - 1] = 3 * (values[n - 1] - values[n - 2]) / h[n - 2] - db;
    } else if (bc == PERIODIC_CONDITION) {
        alpha[0] = 3 * ((values[1] - values[0]) / h[0] - (values[n - 1] - values[n - 2]) / h[n - 2]);
    }

    for (size_t i = 1; i < n - 1; ++i) {
        l[i] = 2 * (intervals[i + 1] - intervals[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    // Step 3: Back substitution to calculate coefficients
    std::vector<double> c(n, 0.0), b(n - 1), d(n - 1);
    for (int j = n - 2; j >= 0; --j) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (values[j + 1] - values[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    // Step 4: Store coefficients in polynomials
    std::vector<Polynomial> polynomials;
    for (size_t i = 0; i < n - 1; ++i) {
        std::vector<double> coefficients = {d[i], c[i], b[i], values[i]};
        polynomials.emplace_back(coefficients);
    }

    piecewisePolynomial = PiecewisePolynomial(polynomials, intervals);
}

PPSpline::PPSpline(const std::vector<double>& intervals, const std::vector<double>& latCoefficients, const std::vector<double>& lonCoefficients, int degree, SplineBoundaryCondition bc, double da, double db)
    : intervals(intervals), degree(degree), da(da), db(db) {
    if (latCoefficients.size() != lonCoefficients.size()) {
        throw std::invalid_argument("Latitude and longitude coefficient sizes must match.");
    }

    // Combine latitude and longitude coefficients into piecewisePolynomial
    for (size_t i = 0; i < latCoefficients.size(); ++i) {
        coefficients.push_back(latCoefficients[i]);
        coefficients.push_back(lonCoefficients[i]);
    }
}
