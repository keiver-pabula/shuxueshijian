#ifndef PP_SPLINE_H
#define PP_SPLINE_H

#include <vector>
#include <ostream>
#include "spline_definitions.h"
#include "../utils/math_function.h"  // Adjust the path to point to the correct location
#include "piecewise_polynomial.h"



class PPSpline {
private:
    std::vector<double> intervals;   // Knot intervals
    std::vector<double> coefficients; // Polynomial coefficients (linear or cubic)
    int degree;                      // Degree of the spline (1 for linear, 3 for cubic)
    double da; // Start derivative
    double db; // End derivative
    PiecewisePolynomial piecewisePolynomial;

    // Helper methods
    void computeLinearSpline(const std::vector<double>& values, int dim);
    void computeCubicSpline(const std::vector<double>& values, int dim, SplineBoundaryCondition bc);
    void computeCubicSplineForPoint2(const std::vector<double>& values, SplineBoundaryCondition bc, double da, double db);    

public:
    // Constructor
    PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree = 1, SplineBoundaryCondition bc = NATURAL_SPLINE);
    PPSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_intervals, SplineBoundaryCondition bc = NATURAL_SPLINE);
    PPSpline(int dim, int order, const std::vector<MathFunction>& functions, double a, double b, int num_intervals, SplineBoundaryCondition bc, bool is_point_2);
    PPSpline(int dim, int order, const std::vector<MathFunction>& functions, const std::vector<double>& time_points, SplineBoundaryCondition bc, double da = 0.0, double db = 0.0);
    PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree, SplineBoundaryCondition bc, double da, double db);
    PPSpline(const std::vector<double>& intervals, const std::vector<double>& latCoefficients, const std::vector<double>& lonCoefficients, int degree, SplineBoundaryCondition bc, double da, double db);
    PPSpline(const std::vector<double>& intervals, const std::vector<double>& values, int degree, SplineBoundaryCondition bc, bool equalSizes);


    // Evaluate the spline at a given point
    double evaluate(double x) const;

    // Print spline details
    void print(std::ostream& os) const;
    void printDetailed(std::ostream &os) const;    

    const std::vector<double>& getCoefficients() const { return coefficients; }
};

#endif // PP_SPLINE_H
