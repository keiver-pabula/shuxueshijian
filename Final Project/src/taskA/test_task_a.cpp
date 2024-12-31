#include "cubic_spline.h"
#include "linear_spline.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip> // For formatted output

double exactFunction(double x) {
    return 1.0 / (1.0 + 25.0 * x * x); // Example: Runge's function
}

std::vector<double> generateNodes(double start, double end, int n) {
    std::vector<double> nodes(n);
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; ++i) {
        nodes[i] = start + i * step;
    }
    return nodes;
}

void exportComparisonData(const std::string& filename, const std::vector<double>& t, const std::vector<double>& splineValues, const std::vector<double>& exactValues) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    outFile << "t,spline,exact\n";
    for (size_t i = 0; i < t.size(); ++i) {
        outFile << t[i] << "," << splineValues[i] << "," << exactValues[i] << "\n";
    }

    outFile.close();
}

double computeMaxError(const std::vector<double>& exactValues, const std::vector<double>& splineValues) {
    double maxError = 0.0;
    for (size_t i = 0; i < exactValues.size(); ++i) {
        double error = std::abs(exactValues[i] - splineValues[i]);
        if (error > maxError) {
            maxError = error;
        }
    }
    return maxError;
}

int main() {
    std::vector<int> subdivisions = {6, 11, 21, 41, 81};
    double start = -1.0, end = 1.0;

    for (int n : subdivisions) {
        std::vector<double> x = generateNodes(start, end, n);
        std::vector<double> y(x.size());
        for (size_t i = 0; i < x.size(); ++i) {
            y[i] = exactFunction(x[i]);
        }

        // Fit cubic spline
        CubicSpline cubicSpline;
        cubicSpline.fitNatural(x, y);

        // Generate finer points for comparison
        std::vector<double> t = generateNodes(start, end, 100);
        std::vector<double> splineValues(t.size()), exactValues(t.size());
        for (size_t i = 0; i < t.size(); ++i) {
            splineValues[i] = cubicSpline.evaluate(t[i]);
            exactValues[i] = exactFunction(t[i]);
        }

        // Export data
        exportComparisonData("comparison_n" + std::to_string(n) + ".csv", t, splineValues, exactValues);

        // Compute and display max error
        double maxError = computeMaxError(exactValues, splineValues);
        std::cout << "Error (N = " << n << "): " << std::setprecision(8) << maxError << std::endl;
    }

    std::cout << "Comparison data exported for visualization." << std::endl;
    return 0;
}
