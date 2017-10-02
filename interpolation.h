//
// Created by ziyinqu on 10/1/17.
//
#include "Eigen/Eigen"
#include "Eigen/src/Core/Matrix.h"

#ifndef MPM_INTERPOLATION_H
#define MPM_INTERPOLATION_H

#endif //MPM_INTERPOLATION_H


Eigen::Vector3f InterpolationKernel(float pos,int& baseNode)
{
    baseNode = std::floor(pos-0.5) + 1;

    Eigen::Vector3f output;

    float d0 = pos - baseNode+1; // 0.5<d0<1.5
    float z = 1.5 - d0;
    output[0] = 0.5 * z * z;

    float d1 = d0 - 1;          // -0.5<d1<0.5
    output[1] = 0.75 - d1 * d1;

    float d2 = 1 - d1;          // 0.5<d2<1.5
    float zz = 1.5 - d2;
    output[2] = 0.5 * zz * zz;

    return output;
}
Eigen::Vector3f GradientInterpolationKernel(float pos,int baseNode)
{
    Eigen::Vector3f graInt;
    float d0 = pos - baseNode + 1;
    float z = 1.5 - d0;

    float d1 = d0 - 1;

    float d2 = 1 - d1;
    float zz = 1.5 - d2;

    graInt[0] = -z;
    graInt[1] = -2*d1;
    graInt[2] = zz;

    return graInt;
}

void QuadraticInterpolation(Eigen::Vector3f particlePointPos, Eigen::Vector3i& baseNode,Eigen::Matrix3f& interpolation,Eigen::Matrix3f& gradientIntp)
{

    Eigen::Vector3f interX = InterpolationKernel(particlePointPos[0],baseNode[0]);
    Eigen::Vector3f interY = InterpolationKernel(particlePointPos[1],baseNode[1]);
    Eigen::Vector3f interZ = InterpolationKernel(particlePointPos[2],baseNode[2]);

    interpolation(0,0)= interX[0];interpolation(0,1) = interX[1];interpolation(0,2)=interX[2];
    interpolation(1,0)= interY[0];interpolation(1,1) = interY[1];interpolation(1,2)=interY[2];
    interpolation(2,0)= interZ[0];interpolation(2,1) = interZ[1];interpolation(2,2)=interZ[2];

    Eigen::Vector3f graIntX = GradientInterpolationKernel(particlePointPos[0],baseNode[0]);
    Eigen::Vector3f graIntY = GradientInterpolationKernel(particlePointPos[1],baseNode[1]);
    Eigen::Vector3f graIntZ = GradientInterpolationKernel(particlePointPos[2],baseNode[2]);

    gradientIntp(0,0)= graIntX[0];gradientIntp(0,1) = graIntX[1];gradientIntp(0,2)=graIntX[2];
    gradientIntp(1,0)= graIntY[0];gradientIntp(1,1) = graIntY[1];gradientIntp(1,2)=graIntY[2];
    gradientIntp(2,0)= graIntZ[0];gradientIntp(2,1) = graIntZ[1];gradientIntp(2,2)=graIntZ[2];
}

