//
// Created by mabc on 2/21/25.
//

#ifndef WAVEPACKET_MODULATION_H
#define WAVEPACKET_MODULATION_H

#include "HelperFunctions.h"
#include <optional>

#include<boost/math/interpolators/cubic_b_spline.hpp>


/**
 * Class to determine the vector potential / electric field that generates the electron modulation.
 * Later this should be extended to solving Maxwells equations near the dielectric?
 * Incoming field is from +x direction.
 */
class EvanescentField {
    double lambda;  // Wavelength of incoming light in nm.
    double I0;  // Maximum intensity of incoming light in 10^9 W/cm^2
    const double c_au = 137.036;  // Speed of light in atomic units
    const double eps0_au = 1./(4.*M_PI);
    const double mu0_au = 1./(eps0_au * std::pow(c_au,2));  // Magnetic permeability in vacuum in atomic units

    // Method for evaluation the derivative of the outgoing Hankel 2 function
    Cdouble Hankel2Deriv(double x, int n) const;

    // Method for calculating coefficients in field expansion
    void calFieldCoeffsPerfectConductor();  // Determine coeffs for perfect conducting cylinder

public:
    double R;  // Radius of cylinder that is scattered upon in nm
    double k;  // Wavevector for the light
    double omega;  //Angular Frequency of incoming light
    double E0, H0;  // Maximum electric and 'magnetic' field amplitude
    cdVec fieldCoeffs;  // Only contains n>=0 since the coeffs are identical for n<0

    EvanescentField(double lambda, double I0, double R);

    // Functions for evaluating the magnetic field - only has a component along z (cylinder axis)
    Cdouble eval_H(double rho, double phi);

    // Functions for evaluating the electric field in different coordinates
    std::pair<Cdouble, Cdouble> eval_E(double rho, double phi) const;  // Get E_rho, E_phi at (rho, phi)
    std::pair<Cdouble, Cdouble> eval_E_dsum(double rho, double phi);  // Brute force evaluate the sum, neglecting symmetry in terms
    std::pair<Cdouble, Cdouble> eval_Exy(double rho, double phi) const;  // Get Ex, Ey at (rho, phi)
    std::pair<Cdouble, Cdouble> eval_CartE(double x, double y) const;  // Get Ex, Ey at (x,y)
    Cdouble eval_CartEx(double x, double y) const;

    // Function for evaluating the asymptotic form of E_rho
    Cdouble asympERho(double rho, double phi);

};

#endif //WAVEPACKET_MODULATION_H
