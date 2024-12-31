#include "b_spline.h"
#include "math_function.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "spline_definitions.h"

// Function to evaluate
double f(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

int main() {
    // Open file streams for cubic and quadratic splines
    std::ofstream cubicFile("output/taskB/cubic.txt");
    std::ofstream quadraticFile("output/taskB/quadratic.txt");

    // Check if files opened successfully
    if (!cubicFile.is_open() || !quadraticFile.is_open()) {
        std::cerr << "Error: Unable to open output files." << std::endl;
        return 1;
    }

    // Define the function to be interpolated
    MathFunction f_func(f);

    // Set up parameters for the domain and interpolation
    double a = -1.0;       // Start of the domain
    double b = 1.0;        // End of the domain
    int num_points = 21;   // Number of intervals in the domain

    // Cubic spline interpolation
    std::vector<MathFunction> f_v_cubic = {f_func};
    BSpline cubicSpline(1, 3, f_v_cubic, a, b, num_points, NATURAL_SPLINE);
    
    cubicFile << "Cubic Spline Interpolation:\n";
    cubicSpline.print(cubicFile);
    cubicFile.close();

    // Quadratic spline interpolation
    std::vector<MathFunction> f_v_quadratic = {f_func};
    BSpline quadraticSpline(1, 2, f_v_quadratic, a, b, num_points, NATURAL_SPLINE);
    
    quadraticFile << "Quadratic Spline Interpolation:\n";
    quadraticSpline.print(quadraticFile);
    quadraticFile.close();

    std::cout << "Cubic and Quadratic spline results written to output/problemB." << std::endl;
    return 0;
}
