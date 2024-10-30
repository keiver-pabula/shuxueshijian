#pragma once
#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <vector>
#include <iostream>

class Interpolation {
public:
    Interpolation(const std::vector<double>& x, const std::vector<double>& y);
    double evaluate(double x) const;

private:
    std::vector<double> x;  
    std::vector<double> y;  
    std::vector<std::vector<double>> dividedDifference; 
};

#endif // INTERPOLATION_H