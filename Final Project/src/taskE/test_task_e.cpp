#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "b_spline.h"

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace TaskESpline {

    // Generate Curve 1: r1(t) = (sqrt(3) * cos(t), 2 / 3 * (sqrt(sqrt(3) * |cos(t)|) + sqrt(3) * sin(t)))
    std::vector<std::vector<double>> generateCurve1(double start, double end, int numPoints) {
        std::vector<std::vector<double>> curve;
        double step = (end - start) / (numPoints - 1);
        for (int i = 0; i < numPoints; ++i) {
            double t = start + i * step;
            curve.push_back({
                sqrt(3) * cos(t),
                2.0 / 3.0 * (sqrt(sqrt(3) * fabs(cos(t))) + sqrt(3) * sin(t))
            });
        }
        return curve;
    }

    // Generate Curve 2: r2(t) = (sin(t) + tcos(t), cos(t) - tsin(t))
    std::vector<std::vector<double>> generateCurve2(double start, double end, int numPoints) {
        std::vector<std::vector<double>> curve;
        double step = (end - start) / (numPoints - 1);
        for (int i = 0; i < numPoints; ++i) {
            double t = start + i * step;
            curve.push_back({
                sin(t) + t * cos(t),
                cos(t) - t * sin(t)
            });
        }
        return curve;
    }

    // Generate Curve 3: r3(t) = (sin(cos(t))cos(sin(t)), sin(cos(t))sin(sin(t)), cos(cos(t)))
    std::vector<std::vector<double>> generateCurve3(double start, double end, int numPoints) {
        std::vector<std::vector<double>> curve;
        double step = (end - start) / (numPoints - 1);
        for (int i = 0; i < numPoints; ++i) {
            double t = start + i * step;
            curve.push_back({
                sin(cos(t)) * cos(sin(t)),
                sin(cos(t)) * sin(sin(t)),
                cos(cos(t))
            });
        }
        return curve;
    }

    // Compute chordal length parameterization

    std::vector<double> computeChordalLength(const std::vector<std::vector<double>>& curve) {
        std::vector<double> params;
        params.push_back(0.0);
        for (size_t i = 1; i < curve.size(); ++i) {
            double dist = 0.0;
            for (size_t j = 0; j < curve[i].size(); ++j) {
                dist += std::pow(curve[i][j] - curve[i - 1][j], 2);
            }
            params.push_back(params.back() + std::sqrt(dist));
        }
        // Normalize to [0, 1]
        double maxParam = params.back();
        for (double& param : params) {
            param /= maxParam;
        }
        return params;
    }

    // Compute uniform parameterization
    std::vector<double> computeUniformParameters(double start, double end, int numPoints) {
        std::vector<double> params;
        double step = (end - start) / (numPoints - 1);
        for (int i = 0; i < numPoints; ++i) {
            params.push_back(start + i * step);
        }
        return params;
    }
}

int main() {
    // Parameters for curves
    std::vector<int> n_values = {10, 40, 160}; // Different resolutions
    std::vector<std::pair<double, double>> curve_ranges = {
        {-M_PI, M_PI},    // Curve 1 range
        {0.0, 6 * M_PI},  // Curve 2 range
        {0.0, 2 * M_PI}   // Curve 3 range
    };

    // Generate curves and fit B-splines for each resolution and parameterization
    for (int resolution : n_values) {
        std::cout << "Processing resolution: N" << resolution << "\n";

        // Curve 1
        auto curve1 = TaskESpline::generateCurve1(curve_ranges[0].first, curve_ranges[0].second, resolution);

        // Chordal parameterization
        auto chordalParams1 = TaskESpline::computeChordalLength(curve1);
        BSpline spline1_chord(std::vector<double>(resolution, 0.0), 3);
        spline1_chord.fit(chordalParams1, chordalParams1); // Fit spline with chordal parameters

        std::ofstream outFile1_chord("output/TaskE/r1_chord_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile1_chord << "x(t) = " << curve1[i][0] << ", y(t) = " << curve1[i][1] << "\n";
        }
        outFile1_chord.close();

        // Uniform parameterization
        auto uniformParams1 = TaskESpline::computeUniformParameters(0.0, 1.0, resolution);
        BSpline spline1_unit(std::vector<double>(resolution, 0.0), 3);
        spline1_unit.fit(uniformParams1, uniformParams1); // Fit spline with uniform parameters

        std::ofstream outFile1_unit("output/TaskE/r1_unit_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile1_unit << "x(t) = " << curve1[i][0] << ", y(t) = " << curve1[i][1] << "\n";
        }
        outFile1_unit.close();

        // Repeat the above steps for Curve 2
        auto curve2 = TaskESpline::generateCurve2(curve_ranges[1].first, curve_ranges[1].second, resolution);

        // Chordal parameterization
        auto chordalParams2 = TaskESpline::computeChordalLength(curve2);
        BSpline spline2_chord(std::vector<double>(resolution, 0.0), 3);
        spline2_chord.fit(chordalParams2, chordalParams2);

        std::ofstream outFile2_chord("output/TaskE/r2_chord_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile2_chord << "x(t) = " << curve2[i][0] << ", y(t) = " << curve2[i][1] << "\n";
        }
        outFile2_chord.close();

        // Uniform parameterization
        auto uniformParams2 = TaskESpline::computeUniformParameters(0.0, 1.0, resolution);
        BSpline spline2_unit(std::vector<double>(resolution, 0.0), 3);
        spline2_unit.fit(uniformParams2, uniformParams2);

        std::ofstream outFile2_unit("output/TaskE/r2_unit_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile2_unit << "x(t) = " << curve2[i][0] << ", y(t) = " << curve2[i][1] << "\n";
        }
        outFile2_unit.close();

        // Repeat the above steps for Curve 3
        auto curve3 = TaskESpline::generateCurve3(curve_ranges[2].first, curve_ranges[2].second, resolution);

        // Chordal parameterization
        auto chordalParams3 = TaskESpline::computeChordalLength(curve3);
        BSpline spline3_chord(std::vector<double>(resolution, 0.0), 3);
        spline3_chord.fit(chordalParams3, chordalParams3);

        std::ofstream outFile3_chord("output/TaskE/r3_chord_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile3_chord << "x(t) = " << curve3[i][0]
                           << ", y(t) = " << curve3[i][1]
                           << ", z(t) = " << curve3[i][2] << "\n";
        }
        outFile3_chord.close();

        // Uniform parameterization
        auto uniformParams3 = TaskESpline::computeUniformParameters(0.0, 1.0, resolution);
        BSpline spline3_unit(std::vector<double>(resolution, 0.0), 3);
        spline3_unit.fit(uniformParams3, uniformParams3);

        std::ofstream outFile3_unit("output/TaskE/r3_unit_N" + std::to_string(resolution) + ".txt");
        for (int i = 0; i < resolution; ++i) {
            outFile3_unit << "x(t) = " << curve3[i][0]
                          << ", y(t) = " << curve3[i][1]
                          << ", z(t) = " << curve3[i][2] << "\n";
        }
        outFile3_unit.close();
    }

    std::cout << "All results have been saved.\n";
    return 0;
}
