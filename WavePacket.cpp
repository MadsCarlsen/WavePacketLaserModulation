//
// Created by au642261 on 3/10/25.
//

#include "WavePacket.h"

double WavePacket::Gauss1D(double k, double k0, double sigma) const {
    const double norm = 1./std::pow((2.*M_PI*sigma*sigma), 1./4.);
    return norm * std::exp(-std::pow(k - k0,2)/(4.*std::pow(sigma,2)));
}

double WavePacket::kzGauss(double kz, double k_perp) const {
    const double perpNorm = 1./std::sqrt((2.*M_PI*sigma_perp*sigma_perp));
    return perpNorm * std::exp(-std::pow(k_perp,2)/(4.*std::pow(sigma_perp,2))) * Gauss1D(kz, k0z, sigma_kz);
}

void WavePacket::setParms_kzGauss(double k0z, double sigma_kz, double sigma_perp) {
    this->k0z = k0z;
    this->sigma_kz = sigma_kz;
    this->sigma_perp = sigma_perp;
}

double WavePacket::EGauss(double k, double theta) const {
    return normEGauss/k * std::exp(-std::pow(k-k0,2)/(4.*std::pow(sigma_k,2))) *
        std::exp(-std::pow(std::sin(theta),2)/(4.*std::pow(sigma_theta,2)));
}

void WavePacket::setParams_EGauss(double k0, double sigma_k, double sigma_theta) {
    this->k0 = k0;
    this->sigma_k = sigma_k;
    this->sigma_theta = sigma_theta;

    // Perform normalization as standard
    normalizeEGauss();
}

void WavePacket::normalizeEGauss() {
    // Normalize the EGauss on the given interval by performing trapz in theta
    double Norm = 2.*M_PI * std::sqrt(2.*M_PI) * sigma_k;  // Phi, k integral

    // Now perform the theta integral
    dVec theta_list = linspace(0, 8*sigma_theta, 2000);  // Some large nr here...
    dVec theta_integrand(theta_list.size());

    for (int i=0; i<theta_list.size(); i++) {
        double sin_i = std::sin(theta_list[i]);
        theta_integrand[i] = sin_i * std::exp(-sin_i*sin_i/(2.*sigma_theta*sigma_theta));
    }
    Norm *= trapz(theta_integrand, theta_list[1]-theta_list[0]);
    normEGauss = 1./std::sqrt(Norm);
}



