#include "cubic_spline.h"
#include "b_spline.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <ostream>
#include <filesystem>

// Exact function for testing (Runge's function)
double exactFunction(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

// Function to generate uniform nodes
std::vector<double> generateNodes(double start, double end, int n) {
    std::vector<double> nodes(n);
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) {
        nodes[i] = start + i * step;
    }
    return nodes;
}

// Export comparison data
void exportComparisonData(const std::string& folder, const std::string& filename, 
                          const std::vector<double>& t,
                          const std::vector<double>& ppFormValues,
                          const std::vector<double>& bSplineValues) {
    std::filesystem::create_directories(folder);
    std::ofstream outFile(folder + "/" + filename);
    if (!outFile) {
        std::cerr << "Error opening file: " << folder + "/" + filename << std::endl;
        return;
    }

    outFile << "t ppForm bSpline\n";
    for (size_t i = 0; i < t.size(); ++i) {
        outFile << t[i] << " " << ppFormValues[i] << " " << bSplineValues[i] << "\n";
    }

    outFile.close();
}

int main() {
    double start = -1.0, end = 1.0;
    int n = 21; // Number of nodes

    // Generate data points
    std::vector<double> x = generateNodes(start, end, n);
    std::vector<double> y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        y[i] = exactFunction(x[i]);
        std::cout << "Node[" << i << "]: x = " << x[i] << ", y = " << y[i] << std::endl;
    }

    // pp-Form (Cubic Spline)
    CubicSpline cubicSpline;
    cubicSpline.fitNatural(x, y);

    // B-Spline (initialize with control points and degree)
    int degree = 3; // Cubic B-spline
    BSpline bSpline(x, degree);
    bSpline.fit(x, y);

    // Get the valid range for evaluation
    double tMin = std::max(cubicSpline.getMinX(), bSpline.getKnots().front());
    double tMax = std::min(cubicSpline.getMaxX(), bSpline.getKnots().back());

    std::cout << "Evaluation Range: [" << tMin << ", " << tMax << "]" << std::endl;

    std::vector<double> t = generateNodes(tMin, tMax, 100);
    std::vector<double> ppFormValues(t.size()), bSplineValues(t.size());

    for (size_t i = 0; i < t.size(); ++i) {
        ppFormValues[i] = cubicSpline.evaluate(t[i]);
        bSplineValues[i] = bSpline.evaluate(t[i]);
        std::cout << "t = " << t[i]
                  << ", pp-Form: " << ppFormValues[i]
                  << ", B-Spline: " << bSplineValues[i] << std::endl;
    }

    // Export comparison data for visualization
    exportComparisonData("output/taskC", "comparison_pp_b_spline.txt", t, ppFormValues, bSplineValues);

    std::cout << "Comparison data exported for Task C.\n";

    return 0;
}
