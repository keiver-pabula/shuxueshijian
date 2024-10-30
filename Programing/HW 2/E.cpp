#include <iostream>
#include <vector>
#include "Interpolation.h"

int main4() {
    std::vector<double> days = { 0, 6, 10, 13, 17, 20, 28 };
    std::vector<double> sample1_weights = { 6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7 };
    std::vector<double> sample2_weights = { 6.67, 16.1, 18.9, 15.0, 10.6, 9.4, 8.89 };

    Interpolation sample1_interpolator(days, sample1_weights);
    Interpolation sample2_interpolator(days, sample2_weights);

    // (a) 
    double future_day = 43;
    double predicted_weight_sample1 = sample1_interpolator.evaluate(future_day);
    double predicted_weight_sample2 = sample2_interpolator.evaluate(future_day);

    std::cout << "Predicted weight on day 43:" << std::endl;
    std::cout << "Sample 1: " << predicted_weight_sample1 << " grams" << std::endl;
    std::cout << "Sample 2: " << predicted_weight_sample2 << " grams" << std::endl;

    // (b) 
    double death_threshold = 5.0; 
    if (predicted_weight_sample1 < death_threshold) {
        std::cout << "Sample 1 is likely to die." << std::endl;
    }
    else {
        std::cout << "Sample 1 is likely to survive." << std::endl;
    }

    if (predicted_weight_sample2 < death_threshold) {
        std::cout << "Sample 2 is likely to die." << std::endl;
    }
    else {
        std::cout << "Sample 2 is likely to survive." << std::endl;
    }

    return 0;
}