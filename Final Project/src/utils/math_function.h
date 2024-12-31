#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include <functional>
#include <stdexcept>

class MathFunction {
private:
    std::function<double(double)> function; // A callable object that represents the function

public:    
    MathFunction() = default;
    
    MathFunction(std::function<double(double)> func) : function(std::move(func)) {
        if (!function) {
            throw std::invalid_argument("Function cannot be null.");
        }
    }

    // Evaluate the function at a specific point
    virtual double evaluate(double x) const {
        if (!function) {
            throw std::runtime_error("Function is not defined.");
        }
        return function(x);
    }

    virtual ~MathFunction() = default;
};

#endif // MATH_FUNCTION_H
