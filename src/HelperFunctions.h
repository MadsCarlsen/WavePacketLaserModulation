//
// Created by au642261 on 2/21/25.
//

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <complex>
#include <iomanip>

using Cdouble = std::complex<double>;
using dVec = std::vector<double>;
using cdVec = std::vector<std::complex<double>>;
using d2Vec = std::vector<std::vector<double>>;
using dPair = std::pair<double, double>;
using cPair = std::pair<Cdouble, Cdouble>;
constexpr std::complex<double> imagI(0., 1.);

// Function template for saving a vector
template <typename T>
void saveVectorToFile(const std::string& filename, const T& vec, int precision=10) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file!\n";
        return;
    }

    outFile << std::setprecision(precision);

    for (auto num : vec) {
        outFile << num << "\n";
    }

    outFile.close();
    std::cout << "Vector saved to " << filename << "\n";
}

// Function template for saving a grid (2D vector)
template <typename T>
void saveGridToFile(const std::string& filename, const T& vec2d) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file!\n";
        return;
    }

    for (auto vec : vec2d) {
        for (auto num : vec) {
            outFile << num << " ";
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Grid saved to " << filename << "\n";
}

// Linspace
dVec linspace(double x_min, double x_max, int Nx);

// Trapz
template<typename T>
T trapz(const std::vector<T>& y, double dx) {
    if (y.size() < 2) return T(0);

    T sum = 0;
    for (size_t i = 0; i < y.size() - 1; ++i) {
        sum += y[i] + y[i+1];
    }
    return 0.5 * dx * sum;
}

// Composite Simpson - Note nr. of intervals must be EVEN, meaning y.size() should be odd
template<typename T>
T simpson(const std::vector<T>& y, double dx) {
    T sum_even = 0;
    T sum_odd = 0;
    int N = y.size();

    for (size_t i = 2; i < N - 1; i += 2) {
        sum_even += y[i];
    }
    for (size_t i = 1; i < N - 1; i += 2) {
        sum_odd += y[i];
    }

    T result = y[0] + y[N - 1] + 2.0 * sum_even + 4.0 * sum_odd;
    return result * dx / 3.0;
}

double trapzInt(const dVec& f_values, double dx);

#endif //HELPERFUNCTIONS_H
