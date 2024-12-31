#ifndef B_SPLINE_H
#define B_SPLINE_H

#include <vector>
#include <string>
#include <ostream>
#include "math_function.h"
#include "spline_definitions.h"
#include "../pp_form/piecewise_polynomial.h"

class BSpline {
private:
    std::vector<double> knots;  // Knot vector
    std::vector<double> controlPoints;  // Control points
    std::vector<PiecewisePolynomial> segments; 
    int degree;  // Degree of the B-spline

    double evaluateBasis(int i, int k, double t) const;  // Evaluate the basis function
    void computeKnots(const std::vector<double>& x);     // Compute the knot vector

public:
    // Constructors
    BSpline(const std::vector<double>& controlPoints, int degree);  // Default constructor
    BSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_points, SplineBoundaryCondition bc);  // Task B constructor
    BSpline(const std::vector<double>& controlPoints, const std::vector<double>& knots, int degree);  // Constructor with explicit knots

    // Getter methods
    std::vector<double> getKnots() const { return knots; }
    int getDegree() const { return degree; }

    // Fit the B-spline to data points
    void fit(const std::vector<double>& x, const std::vector<double>& y);

    // Evaluate the B-spline at a specific t value
    double evaluate(double t) const;

    // Knot computation methods
    void computePeriodicKnots(const std::vector<double>& x); // Periodic boundary condition
    void computeClampedKnots(const std::vector<double>& x);  // Clamped boundary condition
    void computeNaturalKnots(const std::vector<double>& x);  // Natural boundary condition

    // Curve fitting method
    std::vector<double> fitCurve(const std::vector<double>& x, const std::string& method = "uniform");

    // Print method
    void print(std::ostream& os) const;
    void printDetailed(std::ostream& os) const;
    void printDebug(std::ostream& os) const;

};

#endif // B_SPLINE_H
