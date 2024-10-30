#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include "Interpolation.h"

using namespace std;

double exactFunction(double x) {
    return 1.0 / (1 + x * x);
}
vector<double> generateNodes(int n) {
    vector<double> nodes;
    for (int i = 0; i <= n; i++) {
        double x = -5 + 10.0 * i / n;
        nodes.push_back(x);
    }
    return nodes;
}

int main1() {
    vector<int> n_values = { 2, 4, 6, 8 }; 
    ofstream output("interpolation_results.txt");
    if (!output.is_open()) {
        cerr << "无法打开文件进行写入" << endl;
        return 1;
    }

    for (int n : n_values) {
        vector<double> x_values = generateNodes(n);
        vector<double> y_values;

        for (double x : x_values) {
            y_values.push_back(exactFunction(x));
        }

        Interpolation interpolator1(x_values, y_values);

        output << "n = " << n << endl;
        for (double x = -5; x <= 5; x += 0.1) {
            double interpolated = interpolator1.evaluate(x);
            double exact = exactFunction(x);
            output << x << " " << interpolated << " " << exact << endl;
        }
        output << endl;
    }

    output.close();
    cout << "Complete, Saved data to interpolation_results.txt" << endl;
    ifstream input("interpolation_results.txt");
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
