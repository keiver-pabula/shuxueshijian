#pragma once
#include<iostream>
#include<cmath>
#include<limits>
#include "header.h"
using namespace std;

double depth_function(double theta, double R, double h) {
    return R * (1 - cos(theta)) - h;  
}

double depth_function_derivative(double theta, double R) {
    return R * sin(theta);
}

int main5() {
    double R = 10.0; 
    double h = 2.5; 
    double L = 12.4; 
    double tolerance = 1e-6;
    BisectionSolver bisec_solvers;
    NewtonSolver newton_solvers;

    cout << "ʹ�ö��ַ����:" << endl;
    auto depth_eq = [R, h](double theta) { return depth_function(theta, R, h); };
    double theta_bisec = bisec_solvers.solve(depth_eq, 0, 3.14159);
    cout << "���ַ����ý�: " << theta_bisec * (180.0 / 3.14159) << " degrees" << endl;

    cout << "ʹ��ţ�ٷ����:" << endl;
    auto depth_eq_derivative = [R](double theta) { return depth_function_derivative(theta, R); };
    double initial_guess = 1.0; 
    double theta_newton = newton_solvers.solve(depth_eq, initial_guess);
    cout << "ţ�ٷ���ý�: " << theta_newton * (180.0 / 3.14159) << " degrees" << endl;

    return 0;
}