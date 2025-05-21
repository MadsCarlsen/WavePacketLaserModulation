//
// Created by au642261 on 4/25/25.
//

#ifndef WAVEPACKETMODULATION_H
#define WAVEPACKETMODULATION_H

#include "EvanescentField.h"
#include "Grid2D.h"
#include "WavePacket.h"
//#include "Spline2D.h"
#include "HelperFunctions.h"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include "BilinearInterpolatorRegular.h"

using Spline2Dcomplex = BilinearInterpolatorRegular<Cdouble>;

struct ComplexCubicBSpline {
    double xLeft, xRight;
    boost::math::interpolators::cardinal_cubic_b_spline<double> splineReal, splineImag;
    ComplexCubicBSpline(const dVec& xList, const cdVec& funcVals);
    Cdouble eval(double x) const;
};

class WavePacketModulation {
    int maxPhotonOrder;

    public:
    Spline2Dcomplex AkzSpline, AkySpline;
    dVec kzFTlist, kyFTlist;
    Grid2D<Cdouble> AkzGrid, AkyGrid;
    double kzFTmax, kyFTmax, omega, k0, delta_k;
    dVec k_eval;
    WavePacketModulation(double k0, int maxPhotonOrder, EvanescentField& evanescentField, double kzMax, double kyMax, int Nz, int Ny, double impact_param);
    cPair evalAMatElem(double k1, double k2, double theta1, double theta2) const;
    cPair evalAMatElem(double k1, double k2, double theta1, double theta2, const Spline2Dcomplex& Akz, const Spline2Dcomplex& Aky) const;
    Cdouble kIntPV(const cdVec& gList, dVec& kList, double kPole) const;
    std::vector<Grid2D<Cdouble>> evalWPmodCentral(WavePacket psi0, dVec kInitialVec, double kiCentral, dVec theta_i_vec, dVec thetaFinal, int NThetaEval);  // Modulation only calculated for central k, still integrate WP over all ki

    // TODO: Make theta_eval const&?
    void calPhotonAmps(double ki, double theta_i, dVec theta_eval, std::vector<Grid2D<Cdouble>>& alphaGrids, std::vector<Grid2D<Cdouble>>& betaGrids);  // Function for calculating the 'tree' of born amplitudes for given ki, theta_i

    // Calculates the next alpha/beta 'amplitude' functions, assuming that both absorption and emission contributes
    void calNextAlphaBeta(double ki, dVec k_res_vec, dVec theta_eval, Grid2D<Cdouble>& newAlpha, Grid2D<Cdouble>& newBeta,
        const Grid2D<Cdouble>& oldAlpha, const Grid2D<Cdouble>& oldBeta, int li);
    void calFinalAlpha(double ki, double kFinal, const dVec& theta_eval, Grid2D<Cdouble>& newAlpha,
        const Grid2D<Cdouble>& oldAlpha, const Grid2D<Cdouble>& oldBeta, int li);
    void calFinalBeta(double ki, double kFinal, const dVec& theta_eval, Grid2D<Cdouble>& newBeta,
        const Grid2D<Cdouble>& oldAlpha, const Grid2D<Cdouble>& oldBeta, int li);

    // Hacky Bicubic interpolator
    cdVec halfBicubicInterpolate(double k_target, dVec kVec, dVec thetaVec, Grid2D<Cdouble> data);
};


#endif //WAVEPACKETMODULATION_H
