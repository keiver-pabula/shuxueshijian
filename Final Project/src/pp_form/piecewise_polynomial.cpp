#include "piecewise_polynomial.h"
#include <stdexcept>

PiecewisePolynomial::PiecewisePolynomial() = default;



PiecewisePolynomial::PiecewisePolynomial(const std::vector<Polynomial>& polys, const std::vector<double>& pts)
    : polynomials(polys), points(pts) {}

PiecewisePolynomial::PiecewisePolynomial(const std::vector<MathFunction>& latFuncs, 
                                         const std::vector<MathFunction>& lonFuncs, 
                                         const std::vector<double>& points) {
    if (latFuncs.size() != lonFuncs.size() || latFuncs.size() + 1 != points.size()) {
        throw std::invalid_argument("Mismatched sizes for functions and points.");
    }

    for (size_t i = 0; i < latFuncs.size(); ++i) {
        polynomials.push_back(Polynomial({latFuncs[i].evaluate(points[i]), lonFuncs[i].evaluate(points[i])}));
    }
    this->points = points;
}


double PiecewisePolynomial::evaluate(double x) const {
    for (size_t i = 0; i < points.size() - 1; ++i) {
        if (x >= points[i] && x <= points[i + 1]) {
            return polynomials[i].evaluate(x);
        }
    }
    throw std::out_of_range("Value outside range of the piecewise polynomial.");
}

void PiecewisePolynomial::print(std::ostream& os) const {
    for (size_t i = 0; i < polynomials.size(); ++i) {
        // Print the interval
        os << points[i] << "," << points[i + 1] << std::endl;

        // Print coefficients (reverse order to match descending order in output)
        const auto& coeffs = polynomials[i].getCoefficients();
        for (size_t j = coeffs.size(); j > 0; --j) {
            os << coeffs[j - 1] << " ";
        }
        os << std::endl;
    }
}


