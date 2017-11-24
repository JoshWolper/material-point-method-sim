//
// Created by ziyinqu on 11/9/17.
//

#include "interpolation.h"

Eigen::Vector3f calcWeights(float index_space, int& baseNode)
{
    baseNode = std::floor(index_space-0.5);

    Eigen::Vector3f output;

    float d0 = index_space - (float)baseNode; // 0.5<d0<1.5
    float z = 1.5f - d0;
    output[0] = 0.5 * z * z;

    float d1 = d0 - 1;          // -0.5<d1<0.5
    output[1] = 0.75 - d1 * d1;

    float d2 = 1 - d1;          // 0.5<d2<1.5
    float zz = 1.5f - d2;
    output[2] = 0.5 * zz * zz;

    return output;
}

Eigen::Vector4f calcCubicWeights(float index_space, int& baseNode)
{
    baseNode = std::floor(index_space-1);

    Eigen::Vector4f output;

    float d0 = index_space - (float)baseNode; // -2 < d0 < -1
    float z = 2 - d0;
    output[0] =  z * z * z / 6;

    float d1 = d0 - 1;          // -1 < d1 < 0
    output[1] = 0.5 * d1 * d1 * d1 - d1 * d1 + 2/3;

    float d2 = 1 - d1;          //  0 < d2 < 1
    output[2] = 0.5 * d2 * d2 * d2 - d2 * d2 + 2/3;

    float d3 = d2 + 1;          //  1 < d3 < 2
    output[3] = d3 * d3 * d3 / 6;

    return output;
}

Eigen::Vector3f calcGradWeights(float index_space, int baseNode)
{
    Eigen::Vector3f graInt;

    float d0 = index_space - (float)baseNode;
    float z = 1.5f - d0;

    float d1 = d0 - 1;

    float d2 = 1 - d1;
    float zz = 1.5f - d2;

    graInt[0] = -z;
    graInt[1] = -2*d1;
    graInt[2] = zz;

    return graInt;
}

Eigen::Vector4f calcCubicGradWeights(float index_space, int baseNode)
{
    baseNode = std::floor(index_space-1);

    Eigen::Vector4f output;

    float d0 = index_space - (float)baseNode; // -2 < d0 < -1
    float z = 2 - d0;
    output[0] =  0.5 * z * z;

    float d1 = d0 - 1;          // -1 < d1 < 0
    output[1] = 1.5 * d1 * d1 - 2 * d1;

    float d2 = 1 - d1;          //  0 < d2 < 1
    output[2] = 1.5 * d2 * d2 - 2 * d2;

    float d3 = d2 + 1;          //  1 < d3 < 2
    output[3] = 0.5 * d3 * d3;

    return output;

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

void CubicInterpolation(Eigen::Vector3f particlePos, Eigen::Vector3i& baseNode, Eigen::Matrix4f& wp, Eigen::Matrix4f& dwp)
{

    Eigen::Vector4f interX = calcCubicWeights(particlePos[0], baseNode[0]);
    Eigen::Vector4f interY = calcCubicWeights(particlePos[1], baseNode[1]);
    Eigen::Vector4f interZ = calcCubicWeights(particlePos[2], baseNode[2]);

    // calculate weight matrix
    wp(0,0)= interX[0]; wp(0,1) = interX[1]; wp(0,2)=interX[2]; wp(0,3)=interX[3];
    wp(1,0)= interY[0]; wp(1,1) = interY[1]; wp(1,2)=interY[2]; wp(1,3)=interY[3];
    wp(2,0)= interZ[0]; wp(2,1) = interZ[1]; wp(2,2)=interZ[2]; wp(2,3)=interZ[3];

    Eigen::Vector4f graIntX = calcCubicGradWeights(particlePos[0], baseNode[0]);
    Eigen::Vector4f graIntY = calcCubicGradWeights(particlePos[1], baseNode[1]);
    Eigen::Vector4f graIntZ = calcCubicGradWeights(particlePos[2], baseNode[2]);

    // calculate gradient weight matrix
    dwp(0,0)= graIntX[0]; dwp(0,1) = graIntX[1]; dwp(0,2)=graIntX[2]; dwp(0,3)=graIntX[3];
    dwp(1,0)= graIntY[0]; dwp(1,1) = graIntY[1]; dwp(1,2)=graIntY[2]; dwp(1,3)=graIntY[3];
    dwp(2,0)= graIntZ[0]; dwp(2,1) = graIntZ[1]; dwp(2,2)=graIntZ[2]; dwp(2,3)=graIntZ[3];
}