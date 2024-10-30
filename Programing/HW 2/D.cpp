#include <iostream>
#include <vector>
#include "HermiteInterpolation.h"

int main3() {
    std::vector<double> times = { 0, 3, 5, 8, 13 };
    std::vector<double> positions = { 0, 225, 383, 623, 993 };
    std::vector<double> velocities = { 75, 77, 80, 74, 72 };
    HermiteInterpolation interpolator(times, positions, velocities);

    
    double t = 10;
    double predicted_position = interpolator.evaluate(t);
    double predicted_velocity = interpolator.evaluate_diff(t);

    std::cout << "At t = 10 seconds:" << std::endl;
    std::cout << "Predicted Position: " << predicted_position << " feet" << std::endl;
    std::cout << "Predicted Velocity: " << predicted_velocity << " feet/second" << std::endl;

    
    bool exceeded_speed_limit = false;
    for (double t = 0; t <= 13; t += 0.1) {
        double V = interpolator.evaluate_diff(t);
        if (V > 81) {
            std::cout << "The car EXCEEDS the speed limit of 81 feet per second at t = " << t << " seconds." << std::endl;
            exceeded_speed_limit = true;
            break;
        }
    }
    if (!exceeded_speed_limit) {
        std::cout << "The car NEVER exceeds the speed limit of 81 feet per second." << std::endl;
    }

    return 0;
}
