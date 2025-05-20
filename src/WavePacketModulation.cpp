//
// Created by au642261 on 4/25/25.
//

#include "WavePacketModulation.h"
#include "FourierTransformation.h"
#include "Grid2D.h"
#include <algorithm>
#include<boost/math/quadrature/gauss.hpp>
#include <omp.h>

#include <cstdlib>
#include <cmath>


ComplexCubicBSpline::ComplexCubicBSpline(const dVec& xList, const cdVec& funcVals) {
    xLeft = xList[0];
    xRight = xList.back();
    double dx = xList[1] - xLeft;

    // Extract real/imag parts
    dVec func_real, func_imag;
    func_real.reserve(funcVals.size());
    func_imag.reserve(funcVals.size());
    for (Cdouble fi : funcVals) {
        func_real.push_back(fi.real());
        func_imag.push_back(fi.imag());
    }

    // Create the real/imag splines
    splineReal = boost::math::interpolators::cardinal_cubic_b_spline<double> {func_real.begin(), func_real.end(), xLeft, dx};
    splineImag = boost::math::interpolators::cardinal_cubic_b_spline<double> {func_imag.begin(), func_imag.end(), xLeft, dx};
}

Cdouble ComplexCubicBSpline::eval(double x) const{
    // Check if argument is within bounds - if not return 0.
    if (x < xLeft || x > xRight) {throw std::runtime_error("CubicBSpline out of bounds!");};

    // Else return spline values
    return {splineReal(x), splineImag(x)};
}

cdVec WavePacketModulation::halfBicubicInterpolate(double kTarget, dVec kVec, dVec thetaVec, Grid2D<Cdouble> data) {
    // First interpolate along all the kVecs for each theta
    cdVec kSplinedData(thetaVec.size());
    cdVec tempDataSlice(kVec.size());
    for (int i=0; i<thetaVec.size(); i++) {
        // Build the i'th theta spline
        for (int j=0; j<kVec.size(); j++) {
            tempDataSlice[j] = data(j,i);
        }
        // Sample the target k point for each theta
        ComplexCubicBSpline tempSpline(kVec, tempDataSlice);
        kSplinedData[i] = tempSpline.eval(kTarget);
    }
    // Last spline the sampled values over the theta range
    //ComplexCubicBSpline finalSpline(thetaVec, kSplinedData);
    return kSplinedData;
}



WavePacketModulation::WavePacketModulation(double k0, int maxPhotonOrder, EvanescentField& Efield, double kzMax, double kyMax,
                                            int Nz, int Ny, double impact_param) : maxPhotonOrder(maxPhotonOrder), k0(k0) {

    // Grid to evaluate our amplitude functions over?
    omega = Efield.omega;
    delta_k = omega / k0;
    k_eval = linspace(k0 - 4.*delta_k, k0 + 4.*delta_k, 500);  // TODO SHOULD PROBABLY CHECK THE CONVERGENCE OF THIS GUY

    // Create the Fourier transform of the A-fields and save in class as splines
    dVec z_list = FourierTrans::fftRealGrid(Nz, kzMax);
    dVec y_list = FourierTrans::fftRealGrid(Ny, kyMax);
    std::cout << "z_lims FT: " << z_list[0] << " " << z_list.back() << std::endl;

    // Build the spatial A-field in parallel
    std::cout << "Calculating the spatial E-field..." << std::endl;
    Grid2D<Cdouble> Ax_field_spatial(Nz, Ny);
    Grid2D<Cdouble> Ay_field_spatial(Nz, Ny);

    // Interaction limiting Gaussian params
    double sigma_r = 3/Efield.k;
    std::cout << "sigma_r: " << sigma_r << std::endl;

    #pragma omp parallel for
    for (int i=0; i<Nz; i++) {
        for (int j=0; j<Ny; j++) {
            auto [Ex, Ey] = Efield.eval_CartE(-z_list[i], y_list[j] + impact_param);

            // We add Gaussian to limit the interaction region?...
            double r2 = std::pow(z_list[i],2) + std::pow(y_list[j] + impact_param,2);
            double GaussFac = std::exp(-r2/(2*sigma_r*sigma_r));

            Ax_field_spatial(i,j) = -imagI * Ex / Efield.omega * GaussFac;
            Ay_field_spatial(i,j) = -imagI * Ey / Efield.omega * GaussFac;
        }
    }

    // Save these grids?
    //saveVectorToFile("data/z_list.txt", z_list);
    //saveVectorToFile("data/y_list.txt", y_list);
    //Ax_field_spatial.saveToFile("data/Ax.txt");
    //Ay_field_spatial.saveToFile("data/Ay.txt");


    // Prepare the FT
    dVec kz_list = FourierTrans::fftfreq(Nz, z_list[1]-z_list[0]);
    dVec ky_list = FourierTrans::fftfreq(Ny, y_list[1]-y_list[0]);
    kzFTlist = kz_list;
    kyFTlist = ky_list;
    kzFTmax = *std::max_element(kzFTlist.begin(), kzFTlist.end());
    kyFTmax = *std::max_element(kyFTlist.begin(), kyFTlist.end());
    std::cout << kzFTmax << ", " << kyFTmax << std::endl;
    FourierTrans::convert_space_to_FFT(Ax_field_spatial, z_list, y_list, kz_list[0], ky_list[0]);
    FourierTrans::convert_space_to_FFT(Ay_field_spatial, z_list, y_list, kz_list[0], ky_list[0]);

    // Perform the transform
    std::cout << "Performing Fourier transformation to A-field..." << std::endl;
    FourierTrans::fft2d_inplace(Ax_field_spatial);
    FourierTrans::fft2d_inplace(Ay_field_spatial);

    // Convert to momentum function
    double dx = z_list[1] - z_list[0];
    double dy = y_list[1] - y_list[0];
    double dkx = kz_list[1] - kz_list[0];
    double dky = ky_list[1] - ky_list[0];
    FourierTrans::convert_FFT_to_momentum(Ax_field_spatial, z_list[0], y_list[0], dx, dy, dkx, dky);
    FourierTrans::convert_FFT_to_momentum(Ay_field_spatial, z_list[0], y_list[0], dx, dy, dkx, dky);

    // Save these grids?
    //saveVectorToFile("data/kz_list.txt", kz_list);
    //saveVectorToFile("data/ky_list.txt", ky_list);
    //Ax_field_spatial.saveToFile("data/AxFT.txt");
    //Ay_field_spatial.saveToFile("data/AyFT.txt");

    // Make the 2D spline - need to transpose the data first...
    std::cout << "Creating splines of the momentum A-field..." << std::endl;
    AkzGrid = Ax_field_spatial;
    AkyGrid = Ay_field_spatial;
    AkzSpline = Spline2Dcomplex{kz_list, ky_list, Ax_field_spatial};
    AkySpline = Spline2Dcomplex{kz_list, ky_list, Ay_field_spatial};
}

cPair WavePacketModulation::evalAMatElem(double k1, double k2, double theta1, double theta2) const {
    // Eval delta k vector
    double delta_kx = k1*std::cos(theta1) - k2*std::cos(theta2);
    double delta_ky = k1*std::sin(theta1) - k2*std::sin(theta2);

    // Check out of bounds?
    if (std::abs(delta_kx) > kzFTmax || std::abs(delta_ky) > kyFTmax) {
        std::cout << "OUT OF BOUNDS IN AK ELEM EVALUATION!" << std::endl;
        return {imagI*0., imagI*0.};
    }

    // Eval the splines
    Cdouble mel_x = AkzSpline(delta_kx, delta_ky);  // TODO FACTORS OF 2PI!?
    Cdouble mel_y = AkySpline(delta_kx, delta_ky);
    return {mel_x, mel_y};
}

// Same as above but uses local version of splines to enable thread safe iteration?
cPair WavePacketModulation::evalAMatElem(double k1, double k2, double theta1, double theta2, const Spline2Dcomplex& Akz, const Spline2Dcomplex& Aky) const {
    // Eval delta k vector
    double delta_kx = k1*std::cos(theta1) - k2*std::cos(theta2);
    double delta_ky = k1*std::sin(theta1) - k2*std::sin(theta2);

    // Check out of bounds?
    //if (std::abs(delta_kx) > kzFTmax || std::abs(delta_ky) > kyFTmax) {
    //    std::cout << "OUT OF BOUNDS IN AK ELEM EVALUATION!" << std::endl;
    //    return {imagI*0., imagI*0.};
    //}

    // Eval the splines
    Cdouble mel_x = Akz(delta_kx, delta_ky);  // TODO FACTORS OF 2PI!?
    Cdouble mel_y = Aky(delta_kx, delta_ky);
    return {mel_x, mel_y};
}

std::vector<Grid2D<Cdouble>> WavePacketModulation::evalWPmodCentral(WavePacket psi0, dVec kInitialVec, double kiCentral, dVec theta_i_vec, dVec thetaFinal, int NThetaEval) {
    // kInitialVec is the k components of the WP, we only cal. modulation for central ki.
    double delta_theta = 6.* maxPhotonOrder * delta_k / (k0*k0 + k0*delta_k);
    std::vector<Grid2D<Cdouble>> thetaIntegralGrids(2*maxPhotonOrder+1, Grid2D<Cdouble>(theta_i_vec.size(), thetaFinal.size())); // Grid that contains all the theta_i integrands for all theta_f angles and l values of interest

    // Calculate modulation for different angles
    //omp_set_num_threads(10);

    #pragma omp parallel for
    for (int i=0; i<theta_i_vec.size(); i++) {

        if (omp_get_thread_num() == 0) {
            std::cout << "Performing calculation for intial angle: " << theta_i_vec[i]
                   << "    (" << i << "/" << theta_i_vec.size() << ")" << std::endl;
        }

        // For each theta_i we calculate the alpha/beta amplitude functions
        double theta_i = theta_i_vec[i];
        dVec theta_eval = linspace(theta_i - delta_theta, theta_i + delta_theta, NThetaEval);

        // Define alpha/betaGrids, and calculate them for specific ki, theta_i
        std::vector<Grid2D<Cdouble>> alphaGrids;
        std::vector<Grid2D<Cdouble>> betaGrids;
        calPhotonAmps(kiCentral, theta_i, theta_eval, alphaGrids, betaGrids);  // This will be hecking slow
    }
    return thetaIntegralGrids; 
}

        /* 
        // Build each photon order contribution to the WP in the different ki points of the WP
        for (int li=-maxPhotonOrder; li<=maxPhotonOrder; li++) {
            double k_onshell = std::sqrt(std::pow(kiCentral,2) + li*omega);  // Final k for which the WP is evaluated
            std::vector<const Grid2D<Cdouble>*> contributors;  // Will have pointers to the alpha/beta grids that contribute for this l
            cdVec onShellSumVals(theta_eval.size(), imagI*0.);

            // Locate the grids that contribute (except the last order, which is on-shell already)
            for (int ni=std::abs(li); ni<std::abs(maxPhotonOrder); ni++) {
                if ((ni+li)%2 != 0) {continue;}  // These orders have no transition amplitude
                int stateIndex = (ni+li)/2 + (ni-1)*ni/2; // This is beta index, alpha -1

                if (ni == -li) {
                    // Only beta here - add it to lSumGrid
                    contributors.push_back(&betaGrids[stateIndex]);
                } else if (ni == li) {
                    // Only alpha here?
                    contributors.push_back(&alphaGrids[stateIndex-1]);
                } else {
                    // Both alpha and beta contribute here
                    contributors.push_back(&betaGrids[stateIndex]);
                    contributors.push_back(&alphaGrids[stateIndex-1]);
                }
            }
            if (std::abs(li) != maxPhotonOrder) {
                // Sum up all the contributions
                Grid2D<Cdouble> lSumGrid(k_eval.size(),theta_eval.size());
                lSumGrid.fill(Cdouble(0.,0.));  // Initialize to zero
                for (int j=0; j<k_eval.size(); j++) {
                    for (int m=0; m<theta_eval.size(); m++) {
                        for (const Grid2D<Cdouble>* grid : contributors) {
                            lSumGrid(j,m) += (*grid)(j,m);
                        }
                    }
                }
                // Do Bicubic spline to get on-shell value from the summed grids in the theta_eval points
                onShellSumVals = halfBicubicInterpolate(k_onshell, k_eval, theta_eval, lSumGrid);  // TODO CHECK THIS SPLINE IS WORKING OK
            }

            // Last consider the grids for the final order (already on shell)
            int stateIndex = (maxPhotonOrder+li)/2 + (maxPhotonOrder-1)*maxPhotonOrder/2;
            std::vector<const Grid2D<Cdouble>*> finalContributors;
            if ((maxPhotonOrder+li)%2 != 0) {
                // Pass, this Born order has no contribution for this l
            } else if (-li == maxPhotonOrder) {
                // Only beta contribution
                finalContributors.push_back(&betaGrids[stateIndex]);
            } else if (li == maxPhotonOrder) {
                // Only alpha contribution
                finalContributors.push_back(&alphaGrids[stateIndex-1]);
            } else {
                // Both alpha and beta contribution
                finalContributors.push_back(&betaGrids[stateIndex]);
                finalContributors.push_back(&alphaGrids[stateIndex-1]);
            }
            for (int j=0; j<theta_eval.size(); j++) {
                for (const Grid2D<Cdouble>* grid : finalContributors) {
                    onShellSumVals[j] += (*grid)(0,j);
                }
            }

            // Make spline to enable evaluation in thetaFinal points
            ComplexCubicBSpline onShellSpline(theta_eval, onShellSumVals);

            // Evaluate in the thetaFinal points within theta_eval, and save to integrand grid
            for (int j=0; j<thetaFinal.size(); j++) {
                double theta_f = thetaFinal[j];
                if (theta_f > theta_eval[0] & theta_f < theta_eval.back()) {
                    thetaIntegralGrids[li+maxPhotonOrder](i,j) = onShellSpline.eval(theta_f);
                }
            }
            // TODO ADD n=0, l=0 term? What should that be? Just 1?
        }
    }

    // It is time to build our scattered WP!
    std::vector<Grid2D<Cdouble>> finalWPGrid(maxPhotonOrder+1, Grid2D<Cdouble>(kInitialVec.size(), thetaFinal.size()));

    // Get the initial WP on the grid
    const int Nk = static_cast<int>(kInitialVec.size());
    const int Nt = static_cast<int>(theta_i_vec.size());
    Grid2D<double> initialWPGrid(kInitialVec.size(), theta_i_vec.size());
    for (int i=0; i<Nk; i++) {
        for (int j=0; j<Nt; j++) {
            initialWPGrid(i,j) = psi0.EGauss(kInitialVec[i], theta_i_vec[j]);
        }
    }

    // TODO CHECK IF K SHOULD BE IN FINAL INTEGRAL OVER THETA_I!
    // TODO ADD KX DEPENDENCE ALSO (RIGHT NOW PHI=PI/2)
    // Perform integrations for each l wavefunction
    double dTheta = theta_i_vec[1]-theta_i_vec[0];
    for (int li=-maxPhotonOrder; li<=maxPhotonOrder; li++) {
        for (int j=0; j<=kInitialVec.size(); j++) {
            for (int i=0; i<thetaFinal.size(); i++) {
                // Perform the trapz integral on the theta initial
                cdVec theta_i_integrand(theta_i_vec.size());
                for (int m=0; m<theta_i_vec.size(); m++) {
                    theta_i_integrand[i] = initialWPGrid(j,m) * thetaIntegralGrids[li+maxPhotonOrder](m, i);
                }
                finalWPGrid[li+maxPhotonOrder](j,i) = trapz(theta_i_integrand, dTheta);
            }
        }
    }
    return finalWPGrid;
}
*/ 

void WavePacketModulation::calPhotonAmps(double ki, double theta_i, dVec theta_eval, std::vector<Grid2D<Cdouble>>& alphaGrids, std::vector<Grid2D<Cdouble>>& betaGrids) {
    // Build the initial amplitude functions
    Grid2D<Cdouble> alpha_11(k_eval.size(),theta_eval.size());
    Grid2D<Cdouble> beta_11(k_eval.size(),theta_eval.size());
    Grid2D<Cdouble> zeroGrid(k_eval.size(),theta_eval.size());  // Just grid filled with zeros? Kinda hacky.
    zeroGrid.fill(Cdouble(0., 0.));
    for (int m=0; m<k_eval.size(); m++) {
        for (int j=0; j<theta_eval.size(); j++) {
            double km = k_eval[m];
            double theta_j = theta_eval[j];

            // Get the transformed vector potential
            auto [Akz, Aky] = evalAMatElem(km, ki, theta_j, theta_i);
            auto [Akz_c, Aky_c] = evalAMatElem(ki, km, theta_i, theta_j);
            Akz_c = std::conj(Akz_c);
            Aky_c = std::conj(Aky_c);

            // Calculate the dot-product
            alpha_11(m,j) = std::cos(theta_i)*ki*Akz + std::sin(theta_i)*ki*Aky;
            beta_11(m,j) = std::cos(theta_i)*ki*Akz_c + std::sin(theta_i)*ki*Aky_c;
        }
    }
    alphaGrids.push_back(alpha_11);
    betaGrids.push_back(beta_11);

    // Loop over Born orders to build full alpha/beta functions needed for higher orders
    //std::cout << "STARTING CALCULATION OF HIGHER ORDERS!" << std::endl;
    int alphaBetaIndex = 0;
    for (int ni=1; ni<maxPhotonOrder-1; ni++) {
        std::cout << "Current Born order: " << ni << std::endl;
        // Eval next alpha/beta grids using previous ones.
        for (int stateIndex = 0; stateIndex <= ni; stateIndex++) {
            Grid2D<Cdouble> newAlpha(k_eval.size(), theta_eval.size());
            Grid2D<Cdouble> newBeta(k_eval.size(), theta_eval.size());

            std::cout << "    State index: " << stateIndex << std::endl;
            // State index indicates the I_n^l state in a given row in the photon exchange 'triangle'
            int li = -ni + 2*stateIndex;  // Corresponding nr. of exchanged photons

            // Take care of the different possibilities for calculation
            if (stateIndex == 0) {  // Only emission taking place here (only beta func involved)
                calNextAlphaBeta(ki, k_eval, theta_eval, newAlpha, newBeta, zeroGrid, betaGrids[alphaBetaIndex + stateIndex], li);  // <--- GET CORRECT INDEX!
            } else if (stateIndex == ni) {  // Only absorption taking place here (only alpha func involved)
                calNextAlphaBeta(ki, k_eval, theta_eval, newAlpha, newBeta, alphaGrids[alphaBetaIndex + stateIndex], zeroGrid, li);  // <--- GET CORRECT INDEX!
            } else {  // Both - both are good
                calNextAlphaBeta(ki, k_eval, theta_eval, newAlpha, newBeta, alphaGrids[alphaBetaIndex + stateIndex],
                    betaGrids[alphaBetaIndex + stateIndex], li);
            }

            // Save the grids
            alphaGrids.push_back(newAlpha);
            betaGrids.push_back(newBeta);
        }
        alphaBetaIndex += ni;
        // Use the current alpha/beta to add up the on-shell photon exchange amplitudes for all theta_i
    }


    //std::cout << "STARTING CALCULATION OF FINAL ORDER!" << std::endl;
    // Last order we do not need full alpha/beta func, only on-shell k (but still full theta since this is theta_f)
    dVec onShellMomenta(2*maxPhotonOrder+1);
    for (int li=-maxPhotonOrder; li<=maxPhotonOrder; li++) {
        int i = li + maxPhotonOrder;
        onShellMomenta[i] = std::sqrt(std::pow(ki,2) - 2*omega*li);  // TODO PLUS OR MINUS HERE!?
    }

    // Evaluate the final order alpha/beta in the on-shell k points
    for (int stateIndex = 0; stateIndex < maxPhotonOrder; stateIndex++) {
        Grid2D<Cdouble> newAlpha(1, theta_eval.size());
        Grid2D<Cdouble> newBeta(1, theta_eval.size());

        int li = -(maxPhotonOrder-1) + 2*stateIndex;  // Corresponding nr. of photons of the state the new one COMES from
        //std::cout << "    State index: " << stateIndex << std::endl;


        // Take care of the different possibilities for calculation
        // TODO EVALUATE IN THE CORRECT ON SHELL MOMENTA. MEANS THAT WE NEED TO MAKE TWO SEPERATE ALPHA/BETA FUNCTION TO ALLOW DIFFERENT K'S
        if (stateIndex == 0) {  // Only emission taking place here (only beta func involved)
            calFinalBeta(ki, onShellMomenta[li-1+maxPhotonOrder], theta_eval, newBeta, zeroGrid, betaGrids[alphaBetaIndex + stateIndex], li);
            calFinalAlpha(ki, onShellMomenta[li+1+maxPhotonOrder], theta_eval, newAlpha, zeroGrid, betaGrids[alphaBetaIndex + stateIndex], li);
        } else if (stateIndex == maxPhotonOrder-1) {  // Only absorption taking place here (only alpha func involved)
            calFinalBeta(ki, onShellMomenta[li-1+maxPhotonOrder], theta_eval, newBeta, alphaGrids[alphaBetaIndex + stateIndex-1], zeroGrid, li);
            calFinalAlpha(ki, onShellMomenta[li+1+maxPhotonOrder], theta_eval, newAlpha, alphaGrids[alphaBetaIndex + stateIndex-1], zeroGrid, li);  // TODO CHECK THESE INDICES!
        } else {  // Both - both are good
            calFinalBeta(ki, onShellMomenta[li-1+maxPhotonOrder], theta_eval, newBeta, alphaGrids[alphaBetaIndex + stateIndex-1],
                betaGrids[alphaBetaIndex + stateIndex], li);
            calFinalAlpha(ki, onShellMomenta[li+1+maxPhotonOrder], theta_eval, newAlpha, alphaGrids[alphaBetaIndex + stateIndex-1],
                betaGrids[alphaBetaIndex + stateIndex], li);
        }
        // Save the grids
        alphaGrids.push_back(newAlpha);
        betaGrids.push_back(newBeta);
    }
}

Cdouble WavePacketModulation::kIntPV(const cdVec& gList, dVec& kList, double kPole) const {
    double a = kList[0];
    double b = kList.back();
    double dk = kList[1] - kList[0];

    // If pole not on interval - just perform normal trapz
    if (kPole < a || kPole > b) {
        cdVec trapz_list(kList.size());
        for (int i = 0; i < kList.size(); i++) {
            trapz_list[i] = gList[i] / (kList[i] - kPole);
        }
        return trapz<Cdouble>(trapz_list, dk);
    }

    // If on interval - determine index that matches pole and make spline around this index
    int deltaIndex = 20;  // Nr of indices to each side of pole to include in spline?
    int poleIndex = static_cast<int>(std::round((kPole - kList[0])/dk));
    int startIndex = std::clamp(poleIndex - deltaIndex, 0, static_cast<int>(kList.size()) - 1);
    int endIndex   = std::clamp(poleIndex + deltaIndex, 0, static_cast<int>(kList.size()) - 1);

    // Create sublist around the pole
    dVec kSublist(kList.begin() + startIndex, kList.begin() + endIndex + 1);
    cdVec gSublist(gList.begin() + startIndex, gList.begin() + endIndex + 1);

    // Create spline of g-function around the pole
    ComplexCubicBSpline gSpline(kSublist, gSublist);
    Cdouble gPole = gSpline.eval(kPole);

    // Calculate the integrals from start -> spline_start, spline_end -> end
    cdVec intLower(startIndex+1);
    for (int i=0; i<startIndex+1; i++) {
        intLower[i] = (gList[i] - gPole) / (kList[i] - kPole);
    }
    cdVec intUpper(gList.size()-endIndex);
    int j=0;
    for (int i=endIndex; i<gList.size(); i++, j++) {
        intUpper[j] = (gList[i] - gPole) / (kList[i] - kPole);
    }
    Cdouble lower_trapz = trapz(intLower, dk);
    Cdouble upper_trapz = trapz(intUpper, dk);

    /*
    // Save the spline for testing
    cdVec splineTest(kSublist.size());
    for (int i = 0; i < kSublist.size(); i++) {
        splineTest[i] = gSpline.eval(kSublist[i]);
    }
    saveVectorToFile("kSublist.txt", kSublist);
    saveVectorToFile("SubSplineTest.txt", splineTest);
    */

    // Janky way to transition from trapz regular grid, to have pole in center of shifted grid
    int DeltaInt = 4;  // Gauss quade covers DeltaInt grid points to either side - 2*DeltaInt in total
    int NsplineGrid = deltaIndex - DeltaInt;
    dVec kSplineLower = linspace(kList[startIndex], kPole-DeltaInt*dk, NsplineGrid);
    dVec kSplineUpper = linspace(kPole+DeltaInt*dk, kList[endIndex], NsplineGrid);
    cdVec splineLower(NsplineGrid), splineUpper(NsplineGrid);
    for (int i=0; i<NsplineGrid; i++) {
        splineLower[i] = (gSpline.eval(kSplineLower[i]) - gPole) / (kSplineLower[i] - kPole);
        splineUpper[i] = (gSpline.eval(kSplineUpper[i]) - gPole) / (kSplineUpper[i] - kPole);
    }
    Cdouble lowerSplineTrapz = trapz(splineLower, kSplineLower[1]-kSplineLower[0]);
    Cdouble upperSplineTrapz = trapz(splineUpper, kSplineUpper[1]-kSplineUpper[0]);

    // Now determine the integral over the pole by using an even order Gauss-quad
    auto pv_integrand = [DeltaInt, dk, kPole, gSpline, gPole](double x) -> Cdouble {
        return (gSpline.eval(DeltaInt*dk*x + kPole) - gPole) / x;  // Integration interval is substituted to -1 -> 1
    };
    Cdouble central_quad = boost::math::quadrature::gauss<double, 8>::integrate(pv_integrand, -1., 1.);

    // Last include delta-function part of original integral
    Cdouble delta_part = imagI*M_PI*gPole;

    return delta_part + lower_trapz + lowerSplineTrapz + upper_trapz + upperSplineTrapz +
        central_quad + gPole*std::log(std::abs(b-kPole)/std::abs(a-kPole));
}

// Builds next alpha and beta assuming both emission and absorption is happening
void WavePacketModulation::calNextAlphaBeta(double ki, dVec k_res_vec, dVec theta_eval, Grid2D<Cdouble>& newAlpha, Grid2D<Cdouble>& newBeta, const Grid2D<Cdouble> &oldAlpha, const Grid2D<Cdouble> &oldBeta, int li) {
    // li is the photon exchange channel where the old alpha/beta COMES from. New l's are taken care of in the G!
    double dTheta = theta_eval[1] - theta_eval[0];

    // Build the parts of the g-function that does not depend on k_res_vec (everything but the A-term)
    double kl = std::sqrt(std::pow(ki,2) + 2*omega*li);
    Grid2D<Cdouble> g_partial(k_eval.size(), theta_eval.size());
    for (int i=0; i<k_eval.size(); i++) {
        double frac = -2. * k_eval[i] / (k_eval[i] + kl); // This is the negative pole that is never reached
        for (int j = 0; j < theta_eval.size(); j++) {
            g_partial(i, j) = frac * (oldAlpha(i, j) + oldBeta(i, j)); // This is the sum of the old alpha/beta
        }
    }

    // TODO CHECK SOME MORE WHY THERE IS SLOWDOWN FOR MANY CORES...
    //omp_set_num_threads(8);

    #pragma omp parallel
    {
        // Vectors used as storage for integral performed in loops
        cdVec gIntList_abs(k_eval.size());  // The g function for all ks, integrated over theta
        cdVec thetaIntList_abs(theta_eval.size());
        cdVec gIntList_em(k_eval.size());  // The g function for all ks, integrated over theta
        cdVec thetaIntList_em(theta_eval.size());

        // Make spline copy local for each thread? 
        //Grid2D AkzGrid_local = AkzGrid;
        //Grid2D AkyGrid_local = AkyGrid;
        //Grid2D g_partial_local = g_partial;

        //Spline2Dcomplex AkzSpline_local {kzFTlist, kyFTlist, AkzGrid_local};
        //Spline2Dcomplex AkySpline_local {kzFTlist, kyFTlist, AkyGrid_local};

        // Pre-calculate all the trig functions of theta_eval, used to evaluate the matrix elements of A
        dVec sinThetaEval(theta_eval.size());
        dVec cosThetaEval(theta_eval.size());
        for (int i=0; i<theta_eval.size(); i++) {
            sinThetaEval[i] = std::sin(theta_eval[i]);
            cosThetaEval[i] = std::cos(theta_eval[i]);
        }


        // Loop over k and theta for which we want to evaluate the new alpha/beta
        #pragma omp for
        for (int iEval=0; iEval<k_res_vec.size(); iEval++) {
            double k_e = k_res_vec[iEval];

            if (omp_get_thread_num() == 0) {
                std::cout << "iEval: " << iEval << "/" << k_res_vec.size() << std::endl;
            }

            for (int jEval=0; jEval<theta_eval.size(); jEval++) {
                //double theta_e = theta_eval[jEval];

                // Loops to perform integrals
                for (int ip=0; ip<k_eval.size(); ip++) {
                    // Outer loop we perform PV integrals over k
                    double kp = k_eval[ip];

                    for (int jp=0; jp<theta_eval.size(); jp++) {
                        // Inner loop performs theta_integrals
                        //double theta_p = theta_eval[jp];

                        // Evaluate matrix elements of A
                        double delta_kx = k_e * cosThetaEval[jEval] - kp * cosThetaEval[jp];
                        double delta_ky = k_e * sinThetaEval[jEval] - kp * sinThetaEval[jp];
                        Cdouble Akz_abs = AkzSpline(delta_kx, delta_ky);  // TODO FACTORS OF 2PI!?
                        Cdouble Aky_abs = AkySpline(delta_kx, delta_ky);
                        Cdouble Akz_em = AkzSpline(-delta_kx, -delta_ky);
                        Akz_em = std::conj(Akz_em);
                        Cdouble Aky_em = AkySpline(-delta_kx, -delta_ky);
                        Aky_em = std::conj(Aky_em);

                        /*
                        Cdouble Akz_abs = 0.*imagI;
                        Cdouble Aky_abs = 0.*imagI;
                        Cdouble Akz_em = 0.*imagI;
                        Cdouble Aky_em = 0.*imagI;
                        */

                        // Build theta integrand
                        thetaIntList_abs[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_abs + sinThetaEval[jp]*Aky_abs);
                        thetaIntList_em[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_em + sinThetaEval[jp]*Aky_em);
                    }

                    // Integrate up theta result as trapz
                    gIntList_abs[ip] = trapz<Cdouble>(thetaIntList_abs, dTheta);
                    gIntList_em[ip] = trapz<Cdouble>(thetaIntList_em, dTheta);
                }

                // With gThetaList calculated for all kp, perform the principle value integral
                newAlpha(iEval, jEval) = kIntPV(gIntList_abs, k_eval, kl);
                newBeta(iEval, jEval) = kIntPV(gIntList_em, k_eval, kl);

                /*
                // Try saving the gThetaList and then break
                std::cout << "kl: " << kl << std::endl;
                std::cout << "Integration result: " << res(iEval, jEval) << std::endl;
                saveVectorToFile("gThetaTest.txt", gThetaList, 14);
                saveVectorToFile("k_eval.txt", k_eval, 14);
                std::exit(1);
                */
            }
        }
    }
}

// Builds the final alpha function 
void WavePacketModulation::calFinalAlpha(double ki, double kFinal, dVec theta_eval, Grid2D<Cdouble>& newAlpha, const Grid2D<Cdouble> &oldAlpha, const Grid2D<Cdouble> &oldBeta, int li) {
    // li is the photon exchange channel where the old alpha/beta COMES from. New l's are taken care of in the G!
    double dTheta = theta_eval[1] - theta_eval[0];
    size_t NThetaEval = theta_eval.size();


    // Test spline stuff
    double dz = AkzSpline.dx;
    double dy = AkzSpline.dy;
    double dz_inv = AkzSpline.dx_inv;
    double dy_inv = AkzSpline.dy_inv;
    double zMin = AkzSpline.xMin;
    double yMin = AkzSpline.yMin;
    dVec zData = AkzSpline.xData;
    dVec yData = AkzSpline.yData;
    Grid2D<Cdouble> AkzData = AkzSpline.fData;
    Grid2D<Cdouble> AkyData = AkySpline.fData;

    // Build the parts of the g-function that does not depend on k_res_vec (everything but the A-term)
    double kl = std::sqrt(std::pow(ki,2) + 2*omega*li);
    Grid2D<Cdouble> g_partial(k_eval.size(), theta_eval.size());
    for (int i=0; i<k_eval.size(); i++) {
        double frac = -2. * k_eval[i] / (k_eval[i] + kl); // This is the negative pole that is never reached
        for (int j = 0; j < theta_eval.size(); j++) {
            g_partial(i, j) = frac * (oldAlpha(i, j) + oldBeta(i, j)); // This is the sum of the old alpha/beta
        }
    }

    // Vectors used as storage for integral performed in loops
    cdVec gIntList_abs(k_eval.size());  // The g function for all ks, integrated over theta
    cdVec thetaIntList_abs(theta_eval.size());

    // Pre-calculate all the trig functions of theta_eval, used to evaluate the matrix elements of A
    dVec sinThetaEval(theta_eval.size());
    dVec cosThetaEval(theta_eval.size());
    for (int i=0; i<theta_eval.size(); i++) {
        sinThetaEval[i] = std::sin(theta_eval[i]);
        cosThetaEval[i] = std::cos(theta_eval[i]);
    }

    // Loop over k and theta for which we want to evaluate the new alpha/beta
    for (int jEval=0; jEval<theta_eval.size(); jEval++) {

        // Loops to perform integrals
        for (int ip=0; ip<k_eval.size(); ip++) {
            // Outer loop we perform PV integrals over k
            double kp = k_eval[ip];


            // Try to get stuff vectorized
            std::vector<int> ixVec(NThetaEval);
            std::vector<int> iyVec(NThetaEval);
            dVec tVec(NThetaEval);
            dVec uVec(NThetaEval);
            for (int jp=0; jp<NThetaEval; jp++) {
                double delta_kx = kFinal * cosThetaEval[jEval] - kp * cosThetaEval[jp];
                double delta_ky = kFinal * sinThetaEval[jEval] - kp * sinThetaEval[jp];
                ixVec[jp] = std::floor((delta_kx-zMin)*dz_inv);
                iyVec[jp] = std::floor((delta_ky-yMin)*dy_inv);
                tVec[jp] = (delta_kx - (zMin+ixVec[jp]*dz))*dz_inv;  // t,u are between 0,1 on the interval from ix0->ix1.
                uVec[jp] = (delta_ky - (yMin+iyVec[jp]*dy))*dy_inv;
            }

            // Evaluate the spline here - not vectorized?
            for (int jp=0; jp<NThetaEval; jp++) {
                int ix0 = ixVec[jp]; int ix1 = ix0+1;
                int iy0 = iyVec[jp]; int iy1 = iy0+1;
                double t = tVec[jp];
                double u = uVec[jp];
                //double u1 = (1.-u); double t1 = (1.-t);
                //double u1t1 = u1*t1; double tu1 = t*u1; double ut1 = u*t1; double tu = t*u;
                //Cdouble Akz_abs = u1t1*AkzData(ix0,iy0) + tu1*AkzData(ix1,iy0) + tu*AkzData(ix1,iy1) + ut1*AkzData(ix0,iy1);
                //Cdouble Aky_abs = u1t1*AkyData(ix0,iy0) + tu1*AkyData(ix1,iy0) + tu*AkyData(ix1,iy1) + ut1*AkyData(ix0,iy1);
                Cdouble Akz_abs = (1.-t)*(1.-u)*AkzData(ix0,iy0) + t*(1.-u)*AkzData(ix1,iy0) + t*u*AkzData(ix1,iy1) + (1.-t)*u*AkzData(ix0,iy1);
                Cdouble Aky_abs = (1.-t)*(1.-u)*AkyData(ix0,iy0) + t*(1.-u)*AkyData(ix1,iy0) + t*u*AkyData(ix1,iy1) + (1.-t)*u*AkyData(ix0,iy1);
                thetaIntList_abs[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_abs + sinThetaEval[jp]*Aky_abs);
            }


            /*
            for (int jp=0; jp<theta_eval.size(); jp++) {
                // Inner loop performs theta_integrals
                //double theta_p = theta_eval[jp];

                // Evaluate matrix elements of A
                double delta_kx = kFinal * cosThetaEval[jEval] - kp * cosThetaEval[jp];
                double delta_ky = kFinal * sinThetaEval[jEval] - kp * sinThetaEval[jp];
                Cdouble Akz_abs = AkzSpline.evalUnchecked(delta_kx, delta_ky);  // TODO FACTORS OF 2PI!?
                Cdouble Aky_abs = AkySpline.evalUnchecked(delta_kx, delta_ky);

                // Build theta integrand
                thetaIntList_abs[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_abs + sinThetaEval[jp]*Aky_abs);
            }
            */


            // Integrate up theta result as trapz
            gIntList_abs[ip] = trapz<Cdouble>(thetaIntList_abs, dTheta);
        }

        // With gThetaList calculated for all kp, perform the principle value integral
        newAlpha(0, jEval) = kIntPV(gIntList_abs, k_eval, kl);
    }
}

// Builds the final beta function 
void WavePacketModulation::calFinalBeta(double ki, double kFinal, dVec theta_eval, Grid2D<Cdouble>& newBeta, const Grid2D<Cdouble> &oldAlpha, const Grid2D<Cdouble> &oldBeta, int li) {
    // li is the photon exchange channel where the old alpha/beta COMES from. New l's are taken care of in the G!
    double dTheta = theta_eval[1] - theta_eval[0];


    // Test spline stuff
    size_t NThetaEval = theta_eval.size();
    double dz = AkzSpline.dx;
    double dy = AkzSpline.dy;
    double dz_inv = AkzSpline.dx_inv;
    double dy_inv = AkzSpline.dy_inv;
    double zMin = AkzSpline.xMin;
    double yMin = AkzSpline.yMin;
    dVec zData = AkzSpline.xData;
    dVec yData = AkzSpline.yData;
    Grid2D<Cdouble> AkzData = AkzSpline.fData;
    Grid2D<Cdouble> AkyData = AkySpline.fData;


    // Build the parts of the g-function that does not depend on k_res_vec (everything but the A-term)
    double kl = std::sqrt(std::pow(ki,2) + 2*omega*li);
    Grid2D<Cdouble> g_partial(k_eval.size(), theta_eval.size());
    for (int i=0; i<k_eval.size(); i++) {
        double frac = -2. * k_eval[i] / (k_eval[i] + kl); // This is the negative pole that is never reached
        for (int j = 0; j < theta_eval.size(); j++) {
            g_partial(i, j) = frac * (oldAlpha(i, j) + oldBeta(i, j)); // This is the sum of the old alpha/beta
        }
    }

    // Vectors used as storage for integral performed in loops
    cdVec gIntList_em(k_eval.size());  // The g function for all ks, integrated over theta
    cdVec thetaIntList_em(theta_eval.size());

    // Pre-calculate all the trig functions of theta_eval, used to evaluate the matrix elements of A
    dVec sinThetaEval(theta_eval.size());
    dVec cosThetaEval(theta_eval.size());
    for (int i=0; i<theta_eval.size(); i++) {
        sinThetaEval[i] = std::sin(theta_eval[i]);
        cosThetaEval[i] = std::cos(theta_eval[i]);
    }

    // Loop over k and theta for which we want to evaluate the new alpha/beta
    for (int jEval=0; jEval<theta_eval.size(); jEval++) {
        //double theta_e = theta_eval[jEval];

        // Loops to perform integrals
        for (int ip=0; ip<k_eval.size(); ip++) {
            // Outer loop we perform PV integrals over k
            double kp = k_eval[ip];


            // Try to get stuff vectorized
            std::vector<int> ixVec(NThetaEval);
            std::vector<int> iyVec(NThetaEval);
            dVec tVec(NThetaEval);
            dVec uVec(NThetaEval);
            for (int jp=0; jp<NThetaEval; jp++) {
                double delta_kx = -kFinal * cosThetaEval[jEval] + kp * cosThetaEval[jp];
                double delta_ky = -kFinal * sinThetaEval[jEval] + kp * sinThetaEval[jp];
                ixVec[jp] = std::floor((delta_kx-zMin)*dz_inv);
                iyVec[jp] = std::floor((delta_ky-yMin)*dy_inv);
                tVec[jp] = (delta_kx - (zMin+ixVec[jp]*dz))*dz_inv;  // t,u are between 0,1 on the interval from ix0->ix1.
                uVec[jp] = (delta_ky - (yMin+iyVec[jp]*dy))*dy_inv;
            }

            // Evaluate the spline here - not vectorized?
            for (int jp=0; jp<NThetaEval; jp++) {
                int ix0 = ixVec[jp]; int ix1 = ix0+1;
                int iy0 = iyVec[jp]; int iy1 = iy0+1;
                double t = tVec[jp];
                double u = uVec[jp];
                //double u1 = (1.-u); double t1 = (1.-t);
                //double u1t1 = u1*t1; double tu1 = t*u1; double ut1 = u*t1; double tu = t*u;
                //Cdouble Akz_abs = u1t1*AkzData(ix0,iy0) + tu1*AkzData(ix1,iy0) + tu*AkzData(ix1,iy1) + ut1*AkzData(ix0,iy1);
                //Cdouble Aky_abs = u1t1*AkyData(ix0,iy0) + tu1*AkyData(ix1,iy0) + tu*AkyData(ix1,iy1) + ut1*AkyData(ix0,iy1);
                Cdouble Akz_abs = (1.-t)*(1.-u)*AkzData(ix0,iy0) + t*(1.-u)*AkzData(ix1,iy0) + t*u*AkzData(ix1,iy1) + (1.-t)*u*AkzData(ix0,iy1);
                Cdouble Aky_abs = (1.-t)*(1.-u)*AkyData(ix0,iy0) + t*(1.-u)*AkyData(ix1,iy0) + t*u*AkyData(ix1,iy1) + (1.-t)*u*AkyData(ix0,iy1);
                Akz_abs = std::conj(Akz_abs);
                Aky_abs = std::conj(Aky_abs);
                thetaIntList_em[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_abs + sinThetaEval[jp]*Aky_abs);
            }

            /*
            for (int jp=0; jp<theta_eval.size(); jp++) {
                // Inner loop performs theta_integrals
                //double theta_p = theta_eval[jp];

                // Evaluate matrix elements of A
                double delta_kx = kFinal * cosThetaEval[jEval] - kp * cosThetaEval[jp];
                double delta_ky = kFinal * sinThetaEval[jEval] - kp * sinThetaEval[jp];
                Cdouble Akz_em = AkzSpline(-delta_kx, -delta_ky);  // TODO FACTORS OF 2PI!?
                Cdouble Aky_em = AkySpline(-delta_kx, -delta_ky);
                Akz_em = std::conj(Akz_em);
                Aky_em = std::conj(Aky_em);

                // Build theta integrand
                thetaIntList_em[jp] = g_partial(ip, jp) * kp * (cosThetaEval[jp]*Akz_em + sinThetaEval[jp]*Aky_em);
            }
            */

            // Integrate up theta result as trapz
            gIntList_em[ip] = trapz<Cdouble>(thetaIntList_em, dTheta);
        }

        // With gThetaList calculated for all kp, perform the principle value integral
        newBeta(0, jEval) = kIntPV(gIntList_em, k_eval, kl);
    }
}


