#ifndef BEZIER_CURVE_H
#define BEZIER_CURVE_H

#include <vector>
#include <utility>
#include <cmath>

class BezierCurve {
public:
    BezierCurve(const std::vector<std::pair<double, double>>& control_points);
    std::pair<double, double> evaluate(double t) const; 

private:
    std::vector<std::pair<double, double>> control_points;
};


BezierCurve::BezierCurve(const std::vector<std::pair<double, double>>& points)
    : control_points(points) {}

std::pair<double, double> BezierCurve::evaluate(double t) const {
    double x = std::pow(1 - t, 3) * control_points[0].first +
        3 * std::pow(1 - t, 2) * t * control_points[1].first +
        3 * (1 - t) * t * t * control_points[2].first +
        t * t * t * control_points[3].first;

    double y = std::pow(1 - t, 3) * control_points[0].second +
        3 * std::pow(1 - t, 2) * t * control_points[1].second +
        3 * (1 - t) * t * t * control_points[2].second +
        t * t * t * control_points[3].second;

    return { x, y };
}

#endif // BEZIER_CURVE_H