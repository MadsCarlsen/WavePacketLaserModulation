//
// Created by au642261 on 3/10/25.
//

#ifndef WAVEPACKET_H
#define WAVEPACKET_H

#include "HelperFunctions.h"

class WavePacket {
     public:
     double k0z, k0x, k0y;  // Central momenta in z,x,y directions
     double k0;  // Central total momentum
     double sigma_kz, sigma_kx, sigma_ky;  // Momentum widths in z,x,y directions
     double sigma_k;  // Momentum width in total momentum
     double sigma_perp, sigma_theta;  // Perpendicular momentum width and angular width
     double normEGauss=1.;

     // Pointer to WP function
     double (WavePacket::*evalPointer)(double kz, double k_perp) const;

     // Various WPs
     double kzGauss(double kz, double k_perp) const;
     double Gauss1D(double k, double k0, double sigma) const;
     double EGauss(double k, double theta) const;

     // Constructor
     WavePacket(double (WavePacket::*chosenWP)(double kz, double k_perp) const) {evalPointer = chosenWP;}

     // Evaluate function (Pointer magic here to point to correct WP type)
     double evalWP(double kz, double k_perp) const {return (this->*evalPointer)(kz, k_perp);}

     // Functions for setting the various parameters used for different WPs
     void setParms_kzGauss(double k0z, double sigma_kz, double sigma_perp);
     void setParams_EGauss(double k0, double sigma_k, double sigma_theta);
     void normalizeEGauss();  // Calculates and sets the normalization of the E-Gaussian

};

#endif //WAVEPACKET_H
