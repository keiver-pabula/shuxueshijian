#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "header.h"
using namespace std;

double f(double z) {
    return z - tan(z);
}

int main2() {
    double tol = 1e-6;
    int max_iter = 1000;
    NewtonSolver newton_solver;

    double root=newton_solver.solve(f, 4.3);
    std::cout << "f µÄ¸ùÎª: " << root << std::endl;

    return 0;
}