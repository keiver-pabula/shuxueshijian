#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include "header.h"
using namespace std;


double f1(double x) {
    return x * x * x - tan(x);
}

double f2(double x) {
    return x * x - 1;
}

double f3(double x) {
    return x * x + x + 2 * cos(x) - 5;
}

double f4(double x) {
    return x * x * x + 4 * x * x + 3 * x - 5;
}

int main1() {
    double tol = 1e-6;
    int max_iter = 1000;
    BisectionSolver bisec_solver;

    double root1 = bisec_solver.solve(f1, 0, 3);
    double root2 = bisec_solver.solve(f2, 0, 1);
    double root3 = bisec_solver.solve(f3, -3, 2);
    double root4 = bisec_solver.solve(f4, 0, 4);

    std::cout << "f1 的根为: " << root1 << std::endl;
    std::cout << "f2 的根为: " << root2 << std::endl;
    std::cout << "f3 的根为: " << root3 << std::endl;
    std::cout << "f4 的根为: " << root4 << std::endl;

    return 0;
}