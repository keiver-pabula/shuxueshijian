#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <filesystem>
#include "pp_spline.h"
#include "b_spline.h"
#include <cmath>
#include <cstdlib>

// Function to generate test values for `f(x) = e^x - x^2 + 1`
double test_function(double x) {
    return exp(x) - x * x + 1;
}

MathFunction f1(test_function);
std::vector<MathFunction> f_v = {f1};

// Utility function to seed random generator for reproducibility
void seed_random() {
    srand(42); // Fixed seed for reproducibility
}

// Point 1: PP Spline & B-Spline S^0_1
void check_P1() {
    try {
        std::cout << "Checking Point 1: PP Spline & B-Spline S^0_1\n";

        // Create output directories
        std::filesystem::create_directories("Output/Check");
        std::filesystem::create_directories("Figure/Check");

        // Define parameters for pp-spline and b-spline
        PPSpline ppSpline(1, 1, f_v, -1, 1, 40, NATURAL_SPLINE);
        BSpline bspline(1, 1, f_v, -1, 1, 40, NATURAL_SPLINE);

        // Save PP Spline details
        std::ofstream ppFile("Output/Check/P1_ppspline.txt");
        ppSpline.print(ppFile);
        ppFile.close();

        // Save B-Spline details
        std::ofstream bFile("Output/Check/P1_bspline.txt");
        bspline.printDetailed(bFile);
        bFile.close();

    } catch (const std::exception &e) {
        std::cerr << "Error in Point 1: " << e.what() << "\n";
    }
}

// Point 2: PP Spline S^2_3 with 3 different boundary conditions
void check_P2() {
    try {
        std::cout << "Checking Point 2: PP Spline S^2_3 with 3 boundary conditions\n";

        std::vector<double> t1;
        for (int i = 0; i < 11; i++) {
            t1.push_back(-1 + 2.0 * rand() / RAND_MAX); // Random knots in [-1, 1]
        }
        std::sort(t1.begin(), t1.end());

        std::vector<SplineBoundaryCondition> boundary_conditions = {
            NATURAL_SPLINE, CLAMPED, PERIODIC_CONDITION
        };
        std::vector<std::string> bc_names = {
            "natural", "clamped", "periodic"
        };

        for (size_t i = 0; i < boundary_conditions.size(); ++i) {
            PPSpline spline(1, 3, f_v, t1, boundary_conditions[i]);
            std::string filename = "Output/Check/P2_s23_" + bc_names[i] + ".txt";

            // Save spline details
            std::ofstream outFile(filename);
            spline.printDetailed(outFile);
            outFile.close();
        }

    } catch (const std::exception &e) {
        std::cerr << "Error in Point 2: " << e.what() << "\n";
    }
}


// Point 3: B-Spline S^2_3 with 3 different boundary conditions
void check_P3() {
    try {
        std::cout << "Checking Point 3: B-Spline S^2_3 with 3 boundary conditions\n";

        seed_random();

        // Generate random knots in [-1, 1]
        std::vector<double> knots(11);
        for (int i = 0; i < 11; ++i) {
            knots[i] = -1 + 2.0 * rand() / RAND_MAX;
        }
        std::sort(knots.begin(), knots.end());

        std::vector<SplineBoundaryCondition> boundary_conditions = {
            NATURAL_SPLINE, CLAMPED, PERIODIC_CONDITION
        };
        std::vector<std::string> bc_names = {
            "natural", "clamped", "periodic"
        };

        for (size_t i = 0; i < boundary_conditions.size(); ++i) {
            BSpline spline(1, 3, f_v, -1, 1, 40, boundary_conditions[i]);
            std::string filename = "Output/Check/P3_s23_" + bc_names[i] + ".txt";

            // Save spline details
            std::ofstream outFile(filename);
            spline.print(outFile);
            outFile.close();
        }

    } catch (const std::exception &e) {
        std::cerr << "Error in Point 3: " << e.what() << "\n";
    }
}

// Point 4: Compare PP Spline & B-Spline for identical conditions
void check_P4() {
    try {
        std::cout << "Checking Point 4: Compare PP Spline & B-Spline\n";

        // Use the same knots and boundary conditions
        seed_random();

        std::vector<double> knots(11);
        for (int i = 0; i < 11; ++i) {
            knots[i] = -1 + 2.0 * rand() / RAND_MAX;
        }
        std::sort(knots.begin(), knots.end());

        PPSpline ppSpline(1, 3, f_v, knots.front(), knots.back(), knots.size() - 1, NATURAL_SPLINE);
        BSpline bspline(1, 3, f_v, -1, 1, 40, NATURAL_SPLINE);

        // Save spline details
        std::ofstream ppFile("Output/Check/P4_ppspline.txt");
        ppSpline.print(ppFile);
        ppFile.close();

        std::ofstream bFile("Output/Check/P4_bspline.txt");
        bspline.print(bFile);
        bFile.close();

    } catch (const std::exception &e) {
        std::cerr << "Error in Point 4: " << e.what() << "\n";
    }
}

// Point 5: B-Spline in any order
void check_P5() {
    try {
        std::cout << "Checking Point 5: B-Spline in any order\n";

        seed_random();

        // Generate random knots with the correct size
        std::vector<double> knots(14 + 4 + 1);  // Number of control points + degree + 1
        for (size_t i = 0; i < knots.size(); ++i) {
            knots[i] = -0.5 + 2.0 * rand() / RAND_MAX;
        }
        std::sort(knots.begin(), knots.end());

        // Random coefficients
        std::vector<double> coefficients(14);  // Number of control points
        for (size_t i = 0; i < coefficients.size(); ++i) {
            coefficients[i] = -1.0 + 2.0 * rand() / RAND_MAX;
        }

        BSpline spline(coefficients, knots, 4);  // Pass the correct degree (4)

        // Save spline details
        std::ofstream outFile("Output/Check/P5_bspline.txt");
        spline.print(outFile);
        outFile.close();

    } catch (const std::exception &e) {
        std::cerr << "Error in Point 5: " << e.what() << "\n";
    }
}


int main() {
    check_P1();
    check_P2();
    check_P3();
    check_P4();
    check_P5();
    return 0;
}
