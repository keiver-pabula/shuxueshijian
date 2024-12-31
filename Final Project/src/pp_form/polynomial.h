#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <iostream>
#include "../utils/math_function.h"

class Polynomial : public MathFunction {
private:
    std::vector<double> coefficients; // Coefficients of the polynomial (low to high order)

public:
    // Default constructor
    Polynomial();

    // Construct polynomial with coefficients
    Polynomial(const std::vector<double>& coef);

    // Construct polynomial using Newton's interpolation
    Polynomial(const std::vector<double>& x_values, const std::vector<double>& y_values);

    // Overloaded operators for addition, subtraction, multiplication
    Polynomial operator+(const Polynomial& other) const;
    Polynomial operator-(const Polynomial& other) const;
    Polynomial operator*(const Polynomial& other) const;

    // Evaluate the polynomial at point x
    double evaluate(double x) const override;

    // Print the polynomial formula
    void print() const;

    const std::vector<double>& getCoefficients() const;

    // Return the derivative of the polynomial
    Polynomial derivative() const;

    ~Polynomial() override = default;
};

#endif // POLYNOMIAL_H
