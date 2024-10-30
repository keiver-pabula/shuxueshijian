#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include "Interpolation.h"

using namespace std;

vector<double> chebyshevNodes(int n) {
    vector<double> nodes;
    for (int i = 0; i < n; i++) {
        double x = cos((2.0 * i + 1) / (2.0 * n) * 3.14159);
        nodes.push_back(x);
    }
    return nodes;
}

double exactFunction1(double x) {
    return 1.0 / (1 + 25 * x * x);
}

int main2() {
    vector<int> n_values = { 5, 10, 15, 20 }; 

    ofstream output("chebyshev_interpolation_results.txt");
    if (!output.is_open()) {
        cerr << "无法打开文件进行写入" << endl;
        return 1;
    }

    for (int n : n_values) {
        vector<double> x_values = chebyshevNodes(n);
        vector<double> y_values;

        for (double x : x_values) {
            y_values.push_back(exactFunction1(x));
        }

        Interpolation interpolator(x_values, y_values);

        output << "n = " << n << endl;
        for (double x = -1; x <= 1; x += 0.05) {
            double interpolated = interpolator.evaluate(x);
            double exact = exactFunction1(x);
            output << x << " " << interpolated << " " << exact << endl;
        }
        output << endl;
    }

    output.close();
    cout << "Complete, Saved data to  chebyshev_interpolation_results.txt" << endl;

    ifstream input("chebyshev_interpolation_results.txt");
    if (!input.is_open()) {
        cerr << "无法打开文件进行读取" << endl;
        return 1;
    }

    string line;
    while (getline(input, line)) {
        cout << line << endl;
    }
    input.close();



    return 0;
}
