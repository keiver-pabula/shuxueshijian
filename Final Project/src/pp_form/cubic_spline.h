#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <vector>

class CubicSpline {
private:
    std::vector<double> x; // Knots
    std::vector<double> y; // Values at knots
    std::vector<double> a, b, c, d; // Coefficients for each polynomial segment

    void solveTridiagonal(const std::vector<double>& h, const std::vector<double>& alpha);
    
public:
    void fitNatural(const std::vector<double>& x_points, const std::vector<double>& y_points);
    void fitClamped(const std::vector<double>& x_points, const std::vector<double>& y_points, double df0, double dfn);
    void fitPeriodic(const std::vector<double>& x_points, const std::vector<double>& y_points);
    
    void fitCurve(const std::vector<double>& x_points, const std::vector<double>& y_points, bool useChordalLength);

    double evaluate(double t) const;

    double getMinX() const { return x.front(); }
    double getMaxX() const { return x.back(); }
};

#endif // CUBIC_SPLINE_H
