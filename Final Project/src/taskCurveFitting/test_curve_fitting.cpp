#include <filesystem> // For directory creation
#include <fstream>
#include "../curve_fitting/curve_fitting.h"
#include <iomanip>
#include <vector>
#include <cmath> // Include math functions and constants
#include <string>
#include <algorithm> // For std::sort

#ifndef M_PI
#define M_PI 3.14159265358979323846 // Define M_PI if not available
#endif

void saveToTxt(
    const std::string& filename,
    const std::vector<std::vector<double>>& originalPoints,
    const std::vector<std::vector<double>>& planePoints,
    const std::vector<std::vector<double>>& planeSpline,
    const std::vector<std::vector<double>>& sphericalPoints)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    file << std::fixed << std::setprecision(6);

    // Write original points
    file << "original_points: \n";
    for (const auto& point : originalPoints) {
        file << point[0] << "," << point[1] << "," << point[2] << "\n";
    }

    // Write plane points
    file << "plane_points: \n";
    for (const auto& point : planePoints) {
        file << point[0] << "," << point[1] << "," << point[2] << "\n";
    }

    // Write plane spline
    file << "plane_spline: \n";
    for (const auto& point : planeSpline) {
        file << point[0] << "," << point[1] << "," << point[2] << "\n";
    }

    // Write spherical points
    file << "spherical_points: \n";
    for (const auto& point : sphericalPoints) {
        file << point[0] << "," << point[1] << "," << point[2] << "\n";
    }

    file.close();
}

void createOutputDirectories()
{
    std::filesystem::create_directories("output/taskCurveFitting");
    std::filesystem::create_directories("figure/taskCurveFitting");
}

// Function to map spherical coordinates to Cartesian (plane) coordinates
std::vector<double> sphericalToCartesian(const std::vector<double>& spherical) {
    if (spherical.size() != 3) {
        throw std::invalid_argument("Spherical input must have 3 components.");
    }
    double x = 2 * spherical[0] / (2 - spherical[2]);
    double y = 2 * spherical[1] / (2 - spherical[2]);
    return {x, y};
}

// Function to map Cartesian coordinates back to spherical coordinates
std::vector<double> cartesianToSpherical(const std::vector<double>& cartesian) {
    if (cartesian.size() != 2) {
        throw std::invalid_argument("Cartesian input must have 2 components.");
    }
    double x = cartesian[0], y = cartesian[1];
    double k = 4 / (4 + x * x + y * y);
    double x_s = k * x;
    double y_s = k * y;
    double z_s = 2 - 2 * k;
    return {x_s, y_s, z_s};
}

void testCurveFitting()
{
    // Generate random original points on a sphere
    std::vector<std::vector<double>> originalPoints;
    for (int i = 0; i < 10; ++i) {
        double theta = 2 * M_PI * rand() / RAND_MAX; // Azimuthal angle
        double phi = M_PI * rand() / RAND_MAX;      // Polar angle
        double x = sin(phi) * cos(theta);          // X-coordinate
        double y = sin(phi) * sin(theta);          // Y-coordinate
        double z = cos(phi) + 1;                   // Z-coordinate, shifted sphere
        originalPoints.push_back({x, y, z});       // Push {x, y, z} as a single point
    }

    // Convert original points to plane points
    std::vector<std::vector<double>> planePoints;
    for (const auto& point : originalPoints) {
        if (point.size() == 3) { // Ensure the point has 3 components
            auto cartesian = sphericalToCartesian(point);
            planePoints.push_back({cartesian[0], cartesian[1], 0.0});
        }
    }

    // Sort plane points by x-coordinate
    std::sort(planePoints.begin(), planePoints.end(), [](const auto& a, const auto& b) {
        return a[0] < b[0];
    });

    // Extract x and y coordinates for spline fitting
    std::vector<double> x_coords, y_coords;
    for (const auto& point : planePoints) {
        x_coords.push_back(point[0]);
        y_coords.push_back(point[1]);
    }

    // Define boundary derivatives for clamped condition
    double da = 0.0; // First derivative at the start
    double db = 0.0; // First derivative at the end

    // Fit a PPSpline with CLAMPED boundary condition
    PPSpline spline(
        x_coords,    // x-coordinates (intervals)
        y_coords,    // y-coordinates (values)
        3,           // Degree (cubic spline)
        CLAMPED,     // Boundary condition
        true         // Flag to indicate equal sizes
    );

    // Sample the spline for visualization
    std::vector<std::vector<double>> planeSpline;
    double x_start = x_coords.front();
    double x_end = x_coords.back();
    for (int i = 0; i <= 100; ++i) {
        double x = x_start + i * (x_end - x_start) / 100.0;
        double y = spline.evaluate(x);
        planeSpline.push_back({x, y, 0.0});
    }

    // Map plane spline points back to spherical coordinates
    std::vector<std::vector<double>> sphericalPoints;
    for (const auto& point : planeSpline) {
        if (point.size() >= 2) { // Ensure the point has at least 2 components for mapping
            auto spherical = cartesianToSpherical({point[0], point[1]});
            sphericalPoints.push_back(spherical);
        }
    }

    // Save all data to a .txt file
    saveToTxt(
        "output/taskCurveFitting/curve_fitting_output.txt",
        originalPoints,
        planePoints,
        planeSpline,
        sphericalPoints);
}



int main()
{
    createOutputDirectories();

    std::cout << "Testing curve fitting..." << std::endl;
    testCurveFitting();

    return 0;
}
