//
// Created by au642261 on 2/21/25.
//

#include "HelperFunctions.h"

dVec linspace(double x_min, double x_max, int Nx) {
    const double dx = (x_max - x_min) / (Nx-1.);
    dVec x_list {};

    for (int i=0; i<Nx; ++i) {
        x_list.push_back(x_min + dx*i);
    }
    return x_list;
}

double trapzInt(const dVec& f_values, double dx) {
    if (f_values.size() < 2) {
        return 0.0; // Integration requires at least two points
    }
    double integral = 0.0;
    for (size_t i = 0; i < f_values.size() - 1; ++i) {
        integral += (f_values[i] + f_values[i + 1]);
    }
    return integral * 0.5 * dx;
}