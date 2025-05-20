//
// Created by au642261 on 3/20/25.
//

#include "FourierTransformation.h"
#include <fftw3.h>
#include <complex>


void FourierTrans::fft2d_inplace(Grid2D<Cdouble>& grid) {
    int nx = static_cast<int>(grid.nRows());
    int ny = static_cast<int>(grid.nCols());

    // Allocate memory for 2D FFTW complex data (fftw_complex array)
    fftw_complex* fft_data = fftw_alloc_complex(nx * ny);

    // Copy the 2D data into fftw_complex array (flattened)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            fft_data[i * ny + j][0] = grid(i, j).real(); // real part
            fft_data[i * ny + j][1] = grid(i, j).imag(); // imaginary part
        }
    }

    // Create the 2D FFT plan (in-place, forward FFT)
    fftw_plan plan = fftw_plan_dft_2d(
        nx, ny,      // Dimensions of the grid
        fft_data,    // Input (output will be in the same array for in-place)
        fft_data,    // Output
        FFTW_FORWARD,// Forward FFT
        FFTW_ESTIMATE// Planning estimate
    );

    // Execute the FFT
    fftw_execute(plan);

    // Copy back the result into the grid (if necessary)
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            grid(i, j) = std::complex<double>(fft_data[i * ny + j][0], fft_data[i * ny + j][1]);
        }
    }

    // Clean up FFTW plan and memory
    fftw_destroy_plan(plan);
    fftw_free(fft_data);
}

dVec FourierTrans::fftfreq(int N, double dx) {
    double dp = 2.*M_PI/(N*dx);
    double pL = -dp/2. * (N-1.);
    dVec p_list(N);
    for (int i=0; i<N; ++i) {
        p_list[i] = pL + dp*i;
    }
    return p_list;
}

dVec FourierTrans::fftRealGrid(int N, double kMax) {
    double dx = (N-1.)/static_cast<double>(N) * M_PI/kMax;
    double xMax = dx/2. * (N-1.);
    dVec x_list(N);
    for (int i=0; i<N; ++i) {
        x_list[i] = dx*i - xMax;
    }
    return x_list;
}

void FourierTrans::convert_space_to_FFT(Grid2D<Cdouble>& dataGrid, const dVec& x_list, const dVec& y_list, double kx_left, double ky_left) {
    int Nx = static_cast<int>(dataGrid.nRows());
    int Ny = static_cast<int>(dataGrid.nCols());

    // First apply in the x-direction
    for (int i=0; i<Nx; ++i) {
        Cdouble exp_fac = std::exp(-imagI * kx_left * x_list[i]);
        for (int j=0; j<Ny; ++j) {
            dataGrid(i,j) *= exp_fac;
        }
    }

    // Then apply in the y-direction
    for (int i=0; i<Ny; ++i) {
        Cdouble exp_fac = std::exp(-imagI * ky_left * y_list[i]);
        for (int j=0; j<Nx; ++j) {
            dataGrid(j,i) *= exp_fac;
        }
    }
}

void FourierTrans::convert_FFT_to_momentum(Grid2D<Cdouble>& dataGrid, double x_left, double y_left, double dx, double dy, double dkx, double dky) {
    int Nx = static_cast<int>(dataGrid.nRows());
    int Ny = static_cast<int>(dataGrid.nCols());

    // First apply in the x-direction
    for (int i=0; i<Nx; ++i) {
        Cdouble exp_fac = std::exp(-imagI * x_left * dkx * (1.*i));
        for (int j=0; j<Ny; ++j) {
            dataGrid(i,j) *= exp_fac;
        }
    }

    // Then apply in the y-direction
    for (int i=0; i<Ny; ++i) {
        Cdouble exp_fac = std::exp(-imagI * y_left * dky * (1.*i));
        for (int j=0; j<Nx; ++j) {
            dataGrid(j,i) *= exp_fac * dx*dy/(2.*M_PI);  // Include normalization
        }
    }
}

void FourierTrans::applyGaussianWindow(Grid2D<Cdouble> &dataGrid, double sigma_x, double sigma_y=-1.) {
    if (sigma_y < 0.) sigma_y = sigma_x;

    // Get the center indices
    double cx = static_cast<double>(dataGrid.nRows()-1) / 2.0;
    double cy = static_cast<double>(dataGrid.nCols()-1) / 2.0;

    // Apply the Gaussian
    for (int i=0; i<dataGrid.nRows(); ++i) {
        double xi = static_cast<double>(i) - cx;
        for (int j=0; j<dataGrid.nCols(); ++j) {
            double yj = static_cast<double>(j) - cy;
            dataGrid(i,j) *= std::exp(-0.5*std::pow(xi,2)/std::pow(sigma_x,2)
                            -0.5*std::pow(yj,2)/std::pow(sigma_y,2));
        }
    }
}

