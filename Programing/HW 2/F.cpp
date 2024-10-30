#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "BezierCurve.h"

std::vector<std::pair<double, double>> generateHeartPoints(int num_points) {
    std::vector<std::pair<double, double>> points;
    double step = 2.0 * sqrt(3.0) / num_points;

    for (double x = -sqrt(3); x <= sqrt(3); x += step) {
        double y = (2.0 / 3) * (sqrt(3.0 - x * x) + sqrt(fabs(x)));
        points.push_back({ x, y });
        points.push_back({ x, -y }); 
    }
    return points;
}

void approximateHeartCurve(int num_segments, const std::string& filename) {
    std::vector<std::pair<double, double>> heart_points = generateHeartPoints(num_segments);
    std::ofstream output(filename);

    if (!output.is_open()) {
        std::cerr << "无法打开文件进行写入" << std::endl;
        return;
    }

    for (size_t i = 0; i < heart_points.size() - 3; i += 3) {
        std::vector<std::pair<double, double>> control_points = {
            heart_points[i],
            heart_points[i + 1],
            heart_points[i + 2],
            heart_points[i + 3]
        };

        BezierCurve bezier(control_points);

        for (double t = 0; t <= 1; t += 0.01) {
            auto point = bezier.evaluate(t);
            output << point.first << " " << point.second << std::endl;
        }
    }

    output.close();
}

int main() {
    approximateHeartCurve(10, "heart_curve_10.txt");
    approximateHeartCurve(40, "heart_curve_40.txt");
    approximateHeartCurve(160, "heart_curve_160.txt");


    std::cout << "Complete, Saved data to  heart_curve_*.txt 文件中" << std::endl;
    return 0;
}