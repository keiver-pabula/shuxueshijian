#ifndef PIECEWISE_POLYNOMIAL_H
#define PIECEWISE_POLYNOMIAL_H

#include <vector>
#include <iostream>
#include "polynomial.h"

class PiecewisePolynomial {
private:
    std::vector<Polynomial> polynomials; // Polynomials for each segment
    std::vector<double> points;          // Breakpoints (intervals)

public:
    // Default constructor
    PiecewisePolynomial();

    // Constructor with polynomials and points
    PiecewisePolynomial(const std::vector<Polynomial>& polys, const std::vector<double>& pts);
    PiecewisePolynomial(const std::vector<MathFunction>& latFuncs, const std::vector<MathFunction>& lonFuncs, const std::vector<double>& points);

    double evaluate(double x) const;
    void print(std::ostream& os) const;
    size_t size() const { return polynomials.size(); }
    MathFunction getFunction(int i) const {
        if (i < 0 || i >= polynomials.size()) {
            throw std::out_of_range("Invalid function index");
        }
        return MathFunction([=](double x) { return polynomials[i].evaluate(x); });
    };

    const std::vector<double>& getPoints() const { return points; }
    std::vector<double> getValues() const {
        std::vector<double> values;
        for (const auto& poly : polynomials) {
            values.push_back(poly.evaluate(0.0)); // Evaluate at t = 0 for demo purposes
        }
        return values;
    }
    
};

#endif // PIECEWISE_POLYNOMIAL_H
