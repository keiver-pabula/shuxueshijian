#pragma once
#include <iostream>
#include <cmath>
#include <functional>
using namespace std;

// 定义抽象基类
class EquationSolver {
public:
    virtual double solve(std::function<double(double)> f, double a, double b) = 0;
    virtual ~EquationSolver() = default; 
};

// 二分法派生类
class BisectionSolver : public EquationSolver {
private:
    double tolerance = 1e-6; 
    double mid;

public:
    double solve(std::function<double(double)> f, double a, double b) override {
        if (f(a) * f(b) > 0) {
            cout << " No real solution. " << endl;
            return 0.0;
        }

        while ((b - a) >= tolerance) {
            mid = (a + b) / 2.0;

            if (f(mid) == 0.0)
                break;

            else if (f(mid) * f(a) < 0)
                b = mid;
            else
                a = mid;
        }

        cout << "Bisection Method: " << mid << endl; 
        return mid;
    }
};

// 牛顿法派生类
class NewtonSolver : public EquationSolver {
private:
    double tolerance = 1e-6;
    double derivative(std::function<double(double)> f, double x) {
        double h = 1e-6; 
        return (f(x + h) - f(x)) / h;
    }

public:
    double solve(std::function<double(double)> f, double x0, double dummy = 0) override {
        double x = x0;

        while (fabs(f(x)) > tolerance) {
            x = x - f(x) / derivative(f, x); 
        }

        cout << "Newton Method: " << x << endl; 
        return x;
    }
};

// 割线法派生类
class SecantSolver : public EquationSolver {
private:
    double tolerance = 1e-6;
    double x2;

public:
    double solve(std::function<double(double)> f, double x0, double x1) override {
        while (fabs(f(x1)) > tolerance) {
            x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
            x0 = x1;
            x1 = x2;
        }

        cout << "Secant Method: " << x1 << endl;
        return x1;
    }
};