//
// Created by mabc on 2/21/25.
//

#include<boost/math/special_functions/hankel.hpp>
#include<boost/math/special_functions/bessel_prime.hpp>
#include<boost/math/special_functions/bessel.hpp>
#include<boost/math/quadrature/trapezoidal.hpp>
#include<boost/math/quadrature/gauss.hpp>
#include "EvanescentField.h"
#include "FourierTransformation.h"


// Member functions for evanescent field calculations
EvanescentField::EvanescentField(double lambda, double I0, double R)
    : lambda {lambda}, I0 {I0}
{
    // Convert everything into atomic units
    E0 = std::sqrt(I0 *1.0e9 / 3.50945e16);  // Maximum electric field amplitude in a.u.
    omega = 2.*M_PI * c_au / (lambda * 1.0e-9 / 5.29177e-11);
    H0 = c_au/(4.*M_PI) * E0;  // H0 = eps_0 * c * E0
    this->R = 100. * R/5.29177210544;  // Assuming R is in nm

    // Wavevector for the light
    k = omega / c_au;

    // Determine the field expansion constants
    calFieldCoeffsPerfectConductor();
}

Cdouble EvanescentField::Hankel2Deriv(double x, int n) const{
    // Hankel 2 is defined as Jv - i * Yv
    return boost::math::cyl_bessel_j_prime(n, x) - imagI * boost::math::cyl_neumann_prime(n, x);
}

void EvanescentField::calFieldCoeffsPerfectConductor() {
    int max_order = static_cast<int>(k*R) + 10;  // Initial guess on maximal order needed
    fieldCoeffs = {};

    for (int i=0; i<=max_order; i++) {
        double JnDeriv = boost::math::cyl_bessel_j_prime(i, k*R);
        Cdouble Hn2Deriv = Hankel2Deriv(k*R, i);
        fieldCoeffs.push_back(-H0 * JnDeriv / Hn2Deriv);
    }
}

std::pair<Cdouble, Cdouble> EvanescentField::eval_E(double rho, double phi) const{
    if (rho > R) {
        // Outside the cylinder
        double krho = k*rho;

        // Take 0'th order by itself, since no 'negative' term there
        Cdouble Hi = boost::math::cyl_hankel_2(0, krho);
        Cdouble Hip = Hankel2Deriv(krho, 0);
        Cdouble E_rho = 0.*imagI;  // Always 0 for n=0
        Cdouble E_phi = fieldCoeffs[0]*Hip;

        // Sum up the series for E_rho and E_phi (taking account of both +-i terms)
        for (int i=1; i<fieldCoeffs.size(); i++) {
            Hi = boost::math::cyl_hankel_2(i, krho);
            Hip = Hankel2Deriv(krho, i);

            if (i%2 == 0) {
                // i even, E_rho is sin, E_phi is cos
                E_rho += -2.*imagI * fieldCoeffs[i]*Hi * (1.*i) * std::sin(i*(phi-M_PI/2.));
                E_phi += 2. * fieldCoeffs[i]*Hip * std::cos(i*(phi-M_PI/2.));
            } else {
                // i odd, E_rho is cos, E_phi is sin
                E_rho += 2. * (1.*i)  * fieldCoeffs[i]*Hi * std::cos(i*(phi-M_PI/2.));
                E_phi += -2. * imagI * fieldCoeffs[i]*Hip * std::sin(i*(phi-M_PI/2.));
            }
        }

        // Multiply by constants
        E_rho *= -4.*M_PI/(rho*omega);
        E_phi *= 4.*M_PI * imagI * k/omega;

        // Also add incomming wave
        Cdouble E_inc = -k/(omega*eps0_au) * H0 * std::exp(imagI*k*rho*std::cos(phi));
        //E_rho += std::sin(phi) * E_inc;
        //E_phi += std::cos(phi) * E_inc;
        return {E_rho, E_phi};

    } else {
        // Inside the cylinder - this should only run if we are working with dielectric else give 0
        return {imagI*0., imagI*0.};
    }
}

std::pair<Cdouble, Cdouble> EvanescentField::eval_E_dsum(double rho, double phi) {
    double krho = k*rho;

    // Determine the max order from length of coeffs
    int max_order = (fieldCoeffs.size() - 1) / 2;

    // Sum up the series
    Cdouble E_rho = 0.*imagI;
    Cdouble E_phi = 0.*imagI;
    for (int i=-max_order; i<=max_order; i++) {
        int fieldIndex = i + max_order;
        Cdouble Hi = boost::math::cyl_hankel_2(i, krho);
        Cdouble Hip = Hankel2Deriv(krho, i);
        Cdouble exp_fac = std::exp(-imagI*(1.*i)*(phi - M_PI/2.));
        E_rho += (1.*i) * fieldCoeffs[fieldIndex] * exp_fac * Hi;
        E_phi += fieldCoeffs[fieldIndex] * exp_fac * Hip;
    }
    E_rho *= -4.*M_PI/(rho*omega);
    E_phi *= 4.*M_PI * imagI * k/omega;
    return {E_rho, E_phi};

}

std::pair<Cdouble, Cdouble> EvanescentField::eval_Exy(double rho, double phi) const{
    // Get field in cylindrical coordinates
    auto [E_rho, E_phi] = eval_E(rho, phi);

    // Convert to cartesian vector
    Cdouble Ex = E_rho * std::cos(phi) - E_phi * std::sin(phi);
    Cdouble Ey = E_rho * std::sin(phi) + E_phi * std::cos(phi);
    return {Ex, Ey};
}

std::pair<Cdouble, Cdouble> EvanescentField::eval_CartE(double x, double y) const{
    // First convert the cartesian coordinates to cylindrical
    double rho = std::sqrt(x*x + y*y);
    double phi = std::atan2(y, x);
    return eval_Exy(rho, phi);
}

Cdouble EvanescentField::eval_CartEx(double x, double y) const{
    auto [E_x, E_y] = eval_CartE(x, y);
    return E_x;
}

Cdouble EvanescentField::eval_H(double rho, double phi) {
    double krho = k*rho;
    Cdouble Hz = fieldCoeffs[0] * boost::math::cyl_hankel_2(0, krho);

    for (int i=1; i<fieldCoeffs.size(); i++) {
        Cdouble Hi = boost::math::cyl_hankel_2(i, krho);
        if (i%2 == 0) {
            Hz += 2.*fieldCoeffs[i] * Hi * std::cos(i*(phi-M_PI/2.));
        } else {
            Hz += -2.*imagI*fieldCoeffs[i] * Hi * std::sin(i*(phi-M_PI/2.));
        }
    }
    return Hz;
}

Cdouble EvanescentField::asympERho(double rho, double phi) {
    double krho = k*rho;
    Cdouble Esum {imagI*0.};
    for (int i=1; i<fieldCoeffs.size(); i++) {
        Cdouble expfac = std::exp(-imagI*(1.*i)*(phi - M_PI));
        Esum += fieldCoeffs[i] * (-imagI*(1.*i)*expfac + imagI*(1.*i)*std::conj(expfac));
    }
    Cdouble factor = -imagI * std::sqrt(2./(M_PI*k*rho)) / (rho*omega*eps0_au);

    return std::exp(-imagI*(k*rho - M_PI/4.)) * factor;  //* Esum;
}

