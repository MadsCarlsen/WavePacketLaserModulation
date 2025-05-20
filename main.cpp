#include <iostream>
#include "WavePacketModulation.h"
#include "EvanescentField.h"
#include "HelperFunctions.h"
#include "Grid2D.h"
#include "FourierTransformation.h"


int main() {
    // Field and wire settings
    const double R = 100;  // Radius of wire in nm
    const double lambda = 800.;//516.6;  // Wavelength in nm  516.6
    const double I0 = 1.;  // Peak intensity in GW
    EvanescentField Efield{lambda, I0, R};

    // Modulation calculations
    const double k0 = 27.;
    const double omega = Efield.omega;
    const double impact_param = Efield.R + 188.9;  // Approx 10 nm from wire

    // Initial WP params
    const double tau_FWHM = 41.341374575751;  // 1 fs? WP width in time
    double sigma_k = std::sqrt(2.*std::log(2.))/(tau_FWHM*k0);
    const double sigma_theta = 10.e-3;  // Width of WP in angular directions
    const double delta_k = omega / k0;
    const double E0 = std::pow(k0,2)/2.;

    // Usefull prints
    std::cout << "E0: " << std::pow(k0,2)/2. << std::endl;
    std::cout << "Omega: " << omega << std::endl;
    std::cout << "sigma_k: " << sigma_k << std::endl;
    std::cout << "delta_k_p: " << std::sqrt(2.*(E0 + omega)) - std::sqrt(2*E0) << std::endl;
    std::cout << "delta_k_m: " << std::sqrt(2.*(E0 - omega)) - std::sqrt(2*E0) << std::endl;
    std::cout << "delta_k approx: " << omega/k0 << std::endl;

    // Create the initial WP
    WavePacket psi0 {&WavePacket::EGauss};
    psi0.setParams_EGauss(k0, sigma_k, sigma_theta);

    dVec theta_i_vec = linspace(-5.*sigma_theta, 5.*sigma_theta, 10);  // <-- Change this to determine nr of iterations in evalWPmodCentral
    dVec k_i_vec = linspace(k0-5.*sigma_k, k0+5.*sigma_k, 200);
    
    /*
    Grid2D<double> psi0Grid(k_i_vec.size(), theta_i_vec.size());
    for (int i=0; i<k_i_vec.size(); i++) {
        for (int j=0; j<theta_i_vec.size(); j++) {
            psi0Grid(i,j) = psi0.EGauss(k_i_vec[i], theta_i_vec[j]);
        }
    }
    psi0Grid.saveToFile("data/psi0Grid.txt");
    saveVectorToFile("data/theta_i_vec.txt", theta_i_vec);
    saveVectorToFile("data/k_i_vec.txt", k_i_vec);
    */ 

    // Setup for the Fourier transform
    int Nz = 2000;
    int Ny = 2000;
    double kxMax = 0.08;
    double kyMax = 0.08;

    double k_i = k0;
    double theta_i = 0.;//0.5e-5;
    double delta_theta = kxMax/(2.2*k0); //20.* maxPhotonOrder * delta_k / (k0*k0 + k0*delta_k);
    std::cout << "delta_theta: " << delta_theta << std::endl;

    
    // TEST BUILDING THE FULL WP! (saves nothing for now)
    const int maxPhotonOrder = 2;
    WavePacketModulation wpMod{k0, maxPhotonOrder, Efield, kxMax, kyMax, Nz, Ny, impact_param};
    wpMod.k_eval = linspace(k0 - 4. * delta_k, k0 + 4. * delta_k, 500);  // <-- Change this for res in Alpha/Beta calculations
    std::vector<Grid2D<Cdouble>> modWPGrid = wpMod.evalWPmodCentral(psi0, k_i_vec, k0, theta_i_vec, theta_i_vec, 600);


    /*
    // TEST CALCULATION OF HIGHER ORDER (CENTRAL LOOPS ONLY)
    const int maxPhotonOrder = 3;
    WavePacketModulation wpMod{k0, maxPhotonOrder, Efield, kxMax, kyMax, Nz, Ny, impact_param};
    wpMod.k_eval = linspace(k0 - 4. * delta_k, k0 + 4. * delta_k, 500);  // <-- Change this for res in Alpha/Beta calculations
    dVec theta_eval = linspace(theta_i - delta_theta, theta_i + delta_theta, 1000);
    std::vector<Grid2D<Cdouble>> alphaGrids;
    std::vector<Grid2D<Cdouble>> betaGrids;

    wpMod.calPhotonAmps(k_i, theta_i, theta_eval, alphaGrids, betaGrids);

    // For testing just save the alpha/beta grids!
    for (int i=0; i<alphaGrids.size(); i++) {
        std::string filename = "data/alphaGrids_" + std::to_string(i) + ".txt";
        alphaGrids[i].saveToFile(filename);
        filename = "data/betaGrids_" + std::to_string(i) + ".txt";
        betaGrids[i].saveToFile(filename);
    }
    saveVectorToFile("data/theta_eval.txt", theta_eval);
    saveVectorToFile("data/k_eval.txt", wpMod.k_eval);
    */ 




    /* 
    // Eval A-FT on k-grid for tests
    std::cout << "Min kz: " << wpMod.kzFTlist[0] << std::endl;
    std::cout << "Max kz: " << wpMod.kzFTlist.back() << std::endl;
    std::cout << "Min ky: " << wpMod.kyFTlist[0] << std::endl;
    std::cout << "Max ky: " << wpMod.kyFTlist.back() << std::endl;

    dVec kz_eval = linspace(0.001, 0.009, 500);
    dVec ky_eval = linspace(-0.009, 0.009, 500);
    Grid2D<Cdouble> AxFourierRes(kz_eval.size(), ky_eval.size());
    Grid2D<Cdouble> AyFourierRes(kz_eval.size(), ky_eval.size());

    for (int i=0; i<kz_eval.size(); i++) {
        for (int j=0; j<ky_eval.size(); j++) {
            AxFourierRes(i,j) = wpMod.AkzSpline(kz_eval[i], ky_eval[j]);
            AyFourierRes(i,j) = wpMod.AkySpline(kz_eval[i], ky_eval[j]);
        }
    }

    AxFourierRes.saveToFile("AxFourierRes.txt");
    AyFourierRes.saveToFile("AyFourierRes.txt");
    saveVectorToFile("kz_eval.txt", kz_eval);
    saveVectorToFile("ky_eval.txt", ky_eval);
    */ 
}
