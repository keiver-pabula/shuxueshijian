#include <iostream>
#include <cmath>
#include <functional>
#include<limits>
#include "header.h"
using namespace std;

double failure_equation(double alpha, double A, double B, double C, double E) {
    return A * sin(alpha) * cos(alpha) + B * pow(sin(alpha), 2) - C * cos(alpha) - E * sin(alpha);
}

// 方程导数
double failure_derivative(double alpha, double A, double B, double C, double E) {
    return A * (cos(2 * alpha)) + 2 * B * sin(alpha) * cos(alpha) + C * sin(alpha) - E * cos(alpha);
}


int main6() {
    double l = 89;  
    double h = 49;  
    double D1 = 55; 
    double D2 = 30; 
    double beta1 = 11.5 * 3.14159 / 180.0; 
    double tolerance = 1e-6;


    cout << "Part (a) alpha ≈ 33°:" << endl;
    double A1 = l * sin(beta1);
    double B1 = l * cos(beta1);
    double C1 = (h + 0.5 * D1) * sin(beta1) - 0.5 * D1 * tan(beta1);
    double E1 = (h + 0.5 * D1) * cos(beta1) - 0.5 * D1;

    auto f1 = [A1, B1, C1, E1](double alpha) { return failure_equation(alpha, A1, B1, C1, E1); };
    auto df1 = [A1, B1, C1, E1](double alpha) { return failure_derivative(alpha, A1, B1, C1, E1); };
    NewtonSolver newton_solver;

    double alpha_a = newton_solver.solve(f1, 33.0 * 3.1415 / 180.0);
    cout << "Newton Method alpha ≈ " << alpha_a * (180.0 / 3.1415) << " degrees" << endl;


    cout << "Part (b):" << endl;
    double C2 = (h + 0.5 * D2) * sin(beta1) - 0.5 * D2 * tan(beta1);
    double E2 = (h + 0.5 * D2) * cos(beta1) - 0.5 * D2;

    auto f2 = [A1, B1, C2, E2](double alpha) { return failure_equation(alpha, A1, B1, C2, E2); };
    auto df2 = [A1, B1, C2, E2](double alpha) { return failure_derivative(alpha, A1, B1, C2, E2); };

    double alpha_b = newton_solver.solve(f2, 33.0 * 3.1415/ 180.0);
    cout << "Newton Method alpha ≈ " << alpha_b * (180.0 / 3.1415) << " degrees" << endl;

    
    cout << "Part (c):" << endl;
    SecantSolver solver;
    double alpha_c = solver.solve(f2, 10.0 * 3.1415 / 180.0, 60.0 * 3.1415 / 180.0);
    cout << "Secant Method alpha ≈ " << alpha_c * (180.0 / 3.1415) << " degrees" << endl;

    return 0;
}