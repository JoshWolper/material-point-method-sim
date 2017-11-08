//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_INTERPOLATION_H
#define MPM_INTERPOLATION_H

#include "Eigen/Eigen"

Eigen::Vector3f calcWeights(float index_space, int& baseNode)
{
    baseNode = std::floor(index_space-0.5);

    Eigen::Vector3f output;

    float d0 = index_space - baseNode; // 0.5<d0<1.5
    float z = 1.5 - d0;
    output[0] = 0.5 * z * z;

    float d1 = d0 - 1;          // -0.5<d1<0.5
    output[1] = 0.75 - d1 * d1;

    float d2 = 1 - d1;          // 0.5<d2<1.5
    float zz = 1.5 - d2;
    output[2] = 0.5 * zz * zz;

    return output;
}

Eigen::Vector3f calcGradWeights(float index_space, int baseNode)
{
    Eigen::Vector3f graInt;

    float d0 = index_space - baseNode;
    float z = 1.5 - d0;

    float d1 = d0 - 1;

    float d2 = 1 - d1;
    float zz = 1.5 - d2;

    graInt[0] = -z;
    graInt[1] = -2*d1;
    graInt[2] = zz;

    return graInt;
}

void QuadraticInterpolation(Eigen::Vector3f particlePos, Eigen::Vector3i& baseNode, Eigen::Matrix3f& wp, Eigen::Matrix3f& dwp)
{

    Eigen::Vector3f interX = calcWeights(particlePos[0], baseNode[0]);
    Eigen::Vector3f interY = calcWeights(particlePos[1], baseNode[1]);
    Eigen::Vector3f interZ = calcWeights(particlePos[2], baseNode[2]);

    // calculate weight matrix
    wp(0,0)= interX[0]; wp(0,1) = interX[1]; wp(0,2)=interX[2];
    wp(1,0)= interY[0]; wp(1,1) = interY[1]; wp(1,2)=interY[2];
    wp(2,0)= interZ[0]; wp(2,1) = interZ[1]; wp(2,2)=interZ[2];

    Eigen::Vector3f graIntX = calcGradWeights(particlePos[0], baseNode[0]);
    Eigen::Vector3f graIntY = calcGradWeights(particlePos[1], baseNode[1]);
    Eigen::Vector3f graIntZ = calcGradWeights(particlePos[2], baseNode[2]);

    // calculate gradient weight matrix
    dwp(0,0)= graIntX[0]; dwp(0,1) = graIntX[1]; dwp(0,2)=graIntX[2];
    dwp(1,0)= graIntY[0]; dwp(1,1) = graIntY[1]; dwp(1,2)=graIntY[2];
    dwp(2,0)= graIntZ[0]; dwp(2,1) = graIntZ[1]; dwp(2,2)=graIntZ[2];
}

#endif //MPM_INTERPOLATION_H
