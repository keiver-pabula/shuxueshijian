#include "linear_spline.h"
#include <stdexcept>

void LinearSpline::fit(const std::vector<double>& x_points, const std::vector<double>& y_points) {
    if (x_points.size() != y_points.size() || x_points.size() < 2) {
        throw std::invalid_argument("Invalid input sizes for linear spline.");
    }

    x = x_points;
    y = y_points;
}

double LinearSpline::evaluate(double t) const {
    if (t < x.front() || t > x.back()) {
        throw std::out_of_range("Value out of bounds.");
    }

    for (size_t i = 0; i < x.size() - 1; ++i) {
        if (t >= x[i] && t <= x[i + 1]) {
            // Linear interpolation formula
            double slope = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            return y[i] + slope * (t - x[i]);
        }
    }

    throw std::out_of_range("Value out of bounds.");
}
