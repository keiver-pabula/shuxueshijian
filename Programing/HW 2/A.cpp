#include <iostream>
#include <vector>
#include "Interpolation.h"

using namespace std;
Interpolation::Interpolation(const std::vector<double>& x, const std::vector<double>& y) : x(x), y(y) {
    int n = x.size();
    dividedDifference.resize(n, std::vector<double>(n));

    for (int i = 0; i < n; i++) {
        dividedDifference[i][0] = y[i];
    }
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dividedDifference[i][j] = (dividedDifference[i + 1][j - 1] - dividedDifference[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
}

    
double Interpolation::evaluate(double x_value) const {

    double result = dividedDifference[0][0];
    double term = 1.0; 


    for (int i = 1; i < x.size(); i++) {
        term *= (x_value - x[i - 1]);
        result += dividedDifference[0][i] * term;
    }

    return result;
}


