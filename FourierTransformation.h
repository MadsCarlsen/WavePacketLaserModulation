//
// Created by au642261 on 3/20/25.
//

#ifndef FOURIERTRANSFORMATION_H
#define FOURIERTRANSFORMATION_H

#include "Grid2D.h"
#include "HelperFunctions.h"

namespace FourierTrans {
    void fft2d_inplace(Grid2D<Cdouble>& grid);
    void convert_space_to_FFT(Grid2D<Cdouble>& dataGrid, const dVec& x_list, const dVec& y_list, double kx_left, double ky_left);
    void convert_FFT_to_momentum(Grid2D<Cdouble>& dataGrid, double x_left, double y_left, double dx, double dy, double dkx, double dky);
    dVec fftfreq(int Nx, double dx);  // Provides the momentum space grid for a given real space grid
    dVec fftRealGrid(int Nx, double kMax);  // Builds symmetric real space grid with max x such that kMax is max momentum
    void applyGaussianWindow(Grid2D<Cdouble>& dataGrid, double sigma_x, double sigma_y);
}


#endif //FOURIERTRANSFORMATION_H
