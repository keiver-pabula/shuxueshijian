// b_spline_basis.cpp
#include <vector>
#include <iostream>

class BSplineBasis {
public:
    double evaluateBasis(int i, int degree, double t, const std::vector<double>& knots) {
        if (degree == 0) {
            return (t >= knots[i] && t < knots[i + 1]) ? 1.0 : 0.0;
        } else {
            double left = (t - knots[i]) / (knots[i + degree] - knots[i]);
            double right = (knots[i + degree + 1] - t) / (knots[i + degree + 1] - knots[i + 1]);

            double leftTerm = left * evaluateBasis(i, degree - 1, t, knots);
            double rightTerm = right * evaluateBasis(i + 1, degree - 1, t, knots);

            return leftTerm + rightTerm;
        }
    }
};
