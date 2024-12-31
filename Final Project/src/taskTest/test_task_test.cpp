#include "b_spline.h"  // Include your spline implementation
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Define test functions
double f1(double x) {
    return 2 * x + 5;  // Linear function
}

double f2(double x) {
    return x * x + 2 * x + 1;  // Quadratic function
}

double f3(double x) {
    return x * x * x + 2 * x * x + 3 * x + 1;  // Cubic function
}

double f4(double x) {
    return exp(x) * x - log(x);  // Exponential + logarithmic function
}

double f5(double x) {
    return sin(x) * cos(x) / x;  // Trigonometric function
}

// Wrap the functions into MathFunction objects
MathFunction f_func1(f1);
MathFunction f_func2(f2);
MathFunction f_func3(f3);
MathFunction f_func4(f4);
MathFunction f_func5(f5);

int main() {
    // Redirect output to a file
    std::ofstream outFile("output/taskTest/func_results.txt");
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing.\n";
        return 1;
    }

    // Define the list of functions to fit
    std::vector<MathFunction> functionList = {f_func1, f_func2, f_func3, f_func4, f_func5};

    // Define the corresponding ranges for each function
    std::vector<std::pair<double, double>> ranges = {
        {-1.0, 1.0},  // Range for f1
        {-1.0, 1.0},  // Range for f2
        {-1.0, 1.0},  // Range for f3
        {0.1, 1.0},   // Range for f4 (avoid log(0))
        {0.1, 1.0}    // Range for f5 (avoid division by 0)
    };

    // Fit cubic splines for each function and print results
    for (size_t i = 0; i < functionList.size(); ++i) {
        MathFunction f = functionList[i];
        double rangeStart = ranges[i].first;
        double rangeEnd = ranges[i].second;

        outFile << "Fitting cubic spline for function " << i + 1 << " in range [" 
                << rangeStart << ", " << rangeEnd << "]\n";

        // Fit a cubic spline
        BSpline cubicSpline(1, 3, {f}, rangeStart, rangeEnd, 21, NATURAL_SPLINE);
        cubicSpline.print(outFile);

        outFile << "====\n";
    }

    outFile.close();
    std::cout << "Spline fitting results saved to 'output/taskTest/func_results.txt'.\n";

    return 0;
}
