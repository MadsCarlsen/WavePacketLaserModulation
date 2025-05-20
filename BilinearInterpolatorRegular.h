//
// Created by au642261 on 4/30/25.
//

#ifndef BILINEARINTERPOLATORREGULAR_H
#define BILINEARINTERPOLATORREGULAR_H

#include "Grid2D.h"
#include <stdexcept>

#include "HelperFunctions.h"

template <typename T>
class BilinearInterpolatorRegular {
public:
    BilinearInterpolatorRegular() = default;

    BilinearInterpolatorRegular& operator=(const BilinearInterpolatorRegular&) = default;

    BilinearInterpolatorRegular(dVec xData, dVec yData, Grid2D<T> fData)
        : xData(xData), yData(yData), fData(fData)
    {
        xMin = xData[0]; xMax = xData.back();
        yMin = yData[0]; yMax = yData.back();
        dx = xData[1] - xMin; dx_inv = 1./dx;
        dy = yData[1] - yMin; dy_inv = 1./dy;
    }

    // Operator for evaluating the interpolator in a given point
    T operator()(double xPoint, double yPoint) const {
        // Do check to see if point is within limits?
        /*
        if (xPoint < xMin || xPoint >= xMax || yPoint < yMin || yPoint >= yMax) {
            // Raise error
            throw std::out_of_range("Point out of interpolation region in bilinear interpolator!");
        }
        */

        // Determine the grid that surrounds the evaluation point. Assuming REGULAR grid!
        int ix0 = std::floor((xPoint-xMin)*dx_inv);
        int ix1 = ix0 + 1;
        int iy0 = std::floor((yPoint-yMin)*dy_inv);
        int iy1 = iy0 + 1;

        // Eval interpolation
        double t = (xPoint - xData[ix0])*dx_inv;  // t,u are between 0,1 on the interval from ix0->ix1.
        double u = (yPoint - yData[iy0])*dy_inv;
        return (1.-t)*(1.-u)*fData(ix0,iy0) + t*(1.-u)*fData(ix1,iy0) + t*u*fData(ix1,iy1) + (1.-t)*u*fData(ix0,iy1);
    }

    // Same as above but unchecked
    T evalUnchecked(double xPoint, double yPoint) const {
        // Determine the grid that surrounds the evaluation point. Assuming REGULAR grid!
        int ix0 = std::floor((xPoint-xMin)*dx_inv);
        int ix1 = ix0 + 1;
        int iy0 = std::floor((yPoint-yMin)*dy_inv);
        int iy1 = iy0 + 1;

        // Eval interpolation
        double t = (xPoint - xData[ix0])*dx_inv;  // t,u are between 0,1 on the interval from ix0->ix1.
        double u = (yPoint - yData[iy0])*dy_inv;
        return (1.-t)*(1.-u)*fData(ix0,iy0) + t*(1.-u)*fData(ix1,iy0) + t*u*fData(ix1,iy1) + (1.-t)*u*fData(ix0,iy1);
    }


    dVec xData, yData;
    Grid2D<T> fData;
    double xMin, xMax, yMin, yMax;
    double dx, dy, dx_inv, dy_inv;
};

#endif //BILINEARINTERPOLATORREGULAR_H
