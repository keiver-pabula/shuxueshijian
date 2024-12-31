#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <filesystem>
#include "cubic_spline.h"
#include "b_spline.h"

// Define the target function
double f(double x) {
    return 1.0 / (1.0 + x * x);
}

int main() {
    // Create output directories
    std::filesystem::create_directories("output/taskD");
    std::string outputFile = "output/taskD/task_d_spline_comparison.txt";

    // Define nodes for cubic and quadratic splines
    std::vector<double> cubicNodes, quadraticNodes;
    for (int i = -5; i <= 5; ++i) {
        cubicNodes.push_back(i);
    }
    for (int i = -5; i < 5; ++i) {
        quadraticNodes.push_back(i + 0.5);
    }

    // Compute function values
    std::vector<double> cubicValues, quadraticValues;
    for (double node : cubicNodes) {
        cubicValues.push_back(f(node));
    }
    for (double node : quadraticNodes) {
        quadraticValues.push_back(f(node));
    }

    // Fit cubic spline
    CubicSpline cubicSpline;
    cubicSpline.fitNatural(cubicNodes, cubicValues);

    // Fit quadratic B-spline
    BSpline quadraticSpline(quadraticValues, 2);
    quadraticSpline.fit(quadraticNodes, quadraticValues);

    // Determine the valid range for evaluation
    double cubicMin = cubicSpline.getMinX();
    double cubicMax = cubicSpline.getMaxX();
    double quadraticMin = quadraticSpline.getKnots()[quadraticSpline.getDegree()];
    double quadraticMax = quadraticSpline.getKnots()[quadraticSpline.getKnots().size() - quadraticSpline.getDegree() - 1];

    double evalMin = std::max(cubicMin, quadraticMin);
    double evalMax = std::min(cubicMax, quadraticMax);

    // Evaluate splines at uniform points within the valid range
    std::vector<double> evalPoints;
    for (double t = evalMin; t <= evalMax; t += 0.1) {
        evalPoints.push_back(t);
    }

    std::vector<double> cubicResults, quadraticResults;
    for (double t : evalPoints) {
        try {
            cubicResults.push_back(cubicSpline.evaluate(t));
        } catch (const std::exception& e) {
            cubicResults.push_back(NAN); // Use NAN for invalid values
        }

        try {
            quadraticResults.push_back(quadraticSpline.evaluate(t));
        } catch (const std::exception& e) {
            quadraticResults.push_back(NAN); // Use NAN for invalid values
        }
    }

    // Save results to a TXT file for visualization
    std::ofstream outFile(outputFile);
    if (!outFile) {
        std::cerr << "Error: Unable to open file for writing.\n";
        return 1;
    }

    outFile << "t cubic quadratic\n";
    for (size_t i = 0; i < evalPoints.size(); ++i) {
        if (!std::isnan(cubicResults[i]) && !std::isnan(quadraticResults[i])) {
            outFile << evalPoints[i] << " " << cubicResults[i] << " " << quadraticResults[i] << "\n";
        }
    }
    outFile.close();

    std::cout << "Results saved to '" << outputFile << "'.\n";

    // Compute errors at specific points
    std::vector<double> checkPoints = {-3.5, -3, -0.5, 0, 0.5, 3, 3.5};
    for (double point : checkPoints) {
        double cubicError = std::abs(cubicSpline.evaluate(point) - f(point));
        double quadraticError = std::abs(quadraticSpline.evaluate(point) - f(point));
        std::cout << "Error at x = " << point << ":\n";
        std::cout << "  Cubic Spline: " << cubicError << "\n";
        std::cout << "  Quadratic Spline: " << quadraticError << "\n";
    }

    return 0;
}
