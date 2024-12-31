#include "polynomial.h"
#include <cmath>
#include <iomanip>

// Construct polynomial with coefficients
Polynomial::Polynomial(const std::vector<double>& coef) : coefficients(coef) {}

// Construct polynomial using Newton's interpolation
Polynomial::Polynomial(const std::vector<double>& x_values, const std::vector<double>& y_values) {
    int n = x_values.size();
    std::vector<std::vector<double>> divided_diff(n, std::vector<double>(n));
    
    // Fill divided difference table
    for (int i = 0; i < n; ++i) {
        divided_diff[i][0] = y_values[i];
    }
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            divided_diff[i][j] = (divided_diff[i + 1][j - 1] - divided_diff[i][j - 1]) / (x_values[i + j] - x_values[i]);
        }
    }

    // Build polynomial
    coefficients = {divided_diff[0][0]};
    for (int i = 1; i < n; ++i) {
        Polynomial term({1});
        for (int j = 0; j < i; ++j) {
            term = term * Polynomial({-x_values[j], 1});
        }
        *this = *this + term * Polynomial({divided_diff[0][i]});
    }
}

// Overloaded addition operator
Polynomial Polynomial::operator+(const Polynomial& other) const {
    size_t max_size = std::max(coefficients.size(), other.coefficients.size());
    std::vector<double> result(max_size, 0.0);

    for (size_t i = 0; i < coefficients.size(); ++i) {
        result[i] += coefficients[i];
    }
    for (size_t i = 0; i < other.coefficients.size(); ++i) {
        result[i] += other.coefficients[i];
    }

    return Polynomial(result);
}

// Overloaded subtraction operator
Polynomial Polynomial::operator-(const Polynomial& other) const {
    size_t max_size = std::max(coefficients.size(), other.coefficients.size());
    std::vector<double> result(max_size, 0.0);

    for (size_t i = 0; i < coefficients.size(); ++i) {
        result[i] += coefficients[i];
    }
    for (size_t i = 0; i < other.coefficients.size(); ++i) {
        result[i] -= other.coefficients[i];
    }

    return Polynomial(result);
}

// Overloaded multiplication operator
Polynomial Polynomial::operator*(const Polynomial& other) const {
    size_t new_size = coefficients.size() + other.coefficients.size() - 1;
    std::vector<double> result(new_size, 0.0);

    for (size_t i = 0; i < coefficients.size(); ++i) {
        for (size_t j = 0; j < other.coefficients.size(); ++j) {
            result[i + j] += coefficients[i] * other.coefficients[j];
        }
    }

    return Polynomial(result);
}

// Evaluate the polynomial at point x
double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = coefficients.size() - 1; i >= 0; --i) {
        result = result * x + coefficients[i];
    }
    return result;
}

// Return the derivative of the polynomial
Polynomial Polynomial::derivative() const {
    if (coefficients.size() <= 1) {
        return Polynomial({0});
    }

    std::vector<double> result(coefficients.size() - 1);
    for (size_t i = 1; i < coefficients.size(); ++i) {
        result[i - 1] = coefficients[i] * i;
    }

    return Polynomial(result);
}

// Print the polynomial formula
void Polynomial::print() const {
    for (size_t i = coefficients.size(); i > 0; --i) {
        if (coefficients[i - 1] != 0) {
            if (i - 1 > 0) {
                std::cout << std::fixed << std::setprecision(6) << coefficients[i - 1] << "x^" << i - 1 << " ";
            } else {
                std::cout << std::fixed << std::setprecision(6) << coefficients[i - 1];
            }
            if (i > 1 && coefficients[i - 2] >= 0) {
                std::cout << "+ ";
            }
        }
    }
    std::cout << std::endl;
}

const std::vector<double>& Polynomial::getCoefficients() const {
    return coefficients;
}

