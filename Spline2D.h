//
// Created by au642261 on 4/24/25.
//

#ifndef SPLINE2D_H
#define SPLINE2D_H

#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <memory>

class Spline2D {
public:
    Spline2D() = default;

    Spline2D(const std::vector<double>& x_vals,
             const std::vector<double>& y_vals,
             const std::vector<double>& z_vals)
        : x(x_vals), y(y_vals), z(z_vals),
          nx(x_vals.size()), ny(y_vals.size())
    {
        type = gsl_interp2d_bicubic;

        spline.reset(gsl_spline2d_alloc(type, nx, ny));
        xacc.reset(gsl_interp_accel_alloc());
        yacc.reset(gsl_interp_accel_alloc());

        gsl_spline2d_init(spline.get(), x.data(), y.data(), z.data(), nx, ny);
    }

    double eval(double xi, double yi) const {
        if (!spline || !xacc || !yacc)
            throw std::runtime_error("Spline2D not initialized");
        return gsl_spline2d_eval(spline.get(), xi, yi, xacc.get(), yacc.get());
    }

private:
    const gsl_interp2d_type* type = nullptr;

    // Smart pointers with custom deleters
    std::unique_ptr<gsl_spline2d, decltype(&gsl_spline2d_free)> spline{nullptr, gsl_spline2d_free};
    std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> xacc{nullptr, gsl_interp_accel_free};
    std::unique_ptr<gsl_interp_accel, decltype(&gsl_interp_accel_free)> yacc{nullptr, gsl_interp_accel_free};

    std::vector<double> x, y, z;
    size_t nx = 0, ny = 0;
};


// Helper class to handle complex splines. Built on top of the above real 2D spline class
class Spline2Dcomplex {
    Spline2D realSpline, imagSpline;

public:
    Spline2Dcomplex() = default;

    Spline2Dcomplex(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
                    const std::vector<std::complex<double>>& z_vals)
    {
        // 'Unpack' the complex numbers
        std::vector<double> zReal(z_vals.size()), zImag(z_vals.size());
        for (size_t i = 0; i < z_vals.size(); i++) {
            zReal[i] = z_vals[i].real();
            zImag[i] = z_vals[i].imag();
        }

        // Make the two splines
        realSpline = Spline2D(x_vals, y_vals, zReal);
        imagSpline = Spline2D(x_vals, y_vals, zImag);
    };

    std::complex<double> eval(double x, double y) const {
        return {realSpline.eval(x, y), imagSpline.eval(x, y)};
    }
};


#endif //SPLINE2D_H
