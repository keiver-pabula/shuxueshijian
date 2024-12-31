#ifndef LINEAR_SPLINE_H
#define LINEAR_SPLINE_H

#include <vector>

class LinearSpline {
private:
    std::vector<double> x; // Knots
    std::vector<double> y; // Values at knots

public:
    void fit(const std::vector<double>& x_points, const std::vector<double>& y_points);
    double evaluate(double t) const;
};

#endif // LINEAR_SPLINE_H
