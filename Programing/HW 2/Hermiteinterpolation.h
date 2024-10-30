#pragma once
#ifndef HERMITE_INTERPOLATION_H
#define HERMITE_INTERPOLATION_H

#include <vector>
#include <stdexcept>
#include <iostream>

class HermiteInterpolation {
private:
    std::vector<double> x;             
    std::vector<double> y;             
    std::vector<double> dy;           
    std::vector<std::vector<double>> div; 
    const double eps = 1e-7;           

public:
    HermiteInterpolation(const std::vector<double>& xValues, const std::vector<double>& yValues, const std::vector<double>& dyValues) {
        if (xValues.size() != yValues.size() || xValues.size() != dyValues.size()) {
            throw std::invalid_argument("The size of x, y, and dy vectors must be the same.");
        }
        x = xValues;
        y = yValues;
        dy = dyValues;
        int n = 2 * x.size(); 
        div.resize(n, std::vector<double>(n, 0.0));

        for (size_t i = 0; i < x.size(); ++i) {
            div[2 * i][0] = y[i];
            div[2 * i + 1][0] = y[i];
            div[2 * i + 1][1] = dy[i];
            if (i < x.size() - 1) {
                div[2 * i][1] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            }
        }

        for (size_t j = 2; j < n; ++j) {
            for (size_t i = 0; i < n - j; ++i) {
                int idx1 = (i + j) / 2;
                int idx2 = i / 2;
                if (idx1 < x.size() && idx2 < x.size()) {
                    div[i][j] = (div[i + 1][j - 1] - div[i][j - 1]) / (x[idx1] - x[idx2]);
                }
            }
        }
    }


 
    double evaluate(double point) const {
        int n = 2 * x.size();
        double result = div[0][0]; 
        double time = 1.0;      
        for (size_t j = 0; j < n - 1; ++j) {
            time *= (point - x[j / 2]);
            result += div[0][j + 1] * time;
        }
        return result;
    }


    double evaluate_diff(double point) const {
        return (evaluate(point + eps) - evaluate(point - eps)) / (eps*2);
    }
};

#endif // HERMITE_INTERPOLATION_H
