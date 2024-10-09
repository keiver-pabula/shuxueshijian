#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "header.h"
using namespace std;

double f5(double x) {
    return sin(x / 2) - 1;
}
double f6(double x) {
    return exp(x) - tan(x);
}
double f7(double x) {
    return pow(x, 3) - 12 * pow(x, 2) + 3 * x + 1;
}

int main() {
    SecantSolver solver;
    double pi = 3.141592653589793;
    solver.solve(f5, 0, pi / 2);
    solver.solve(f5, 1, pi / 2);
    solver.solve(f6, 1, 1.4);
    solver.solve(f6, 0, 1.4);
    solver.solve(f7, 0, -0.5);
    solver.solve(f7, 1, -0.5);
    return 0;
}