//
// Created by VElysianP on 11/7/2017.
//

#ifndef MPM_SVD_H
#define MPM_SVD_H

#endif //MPM_SVD_H

#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/src/SVD/JacobiSVD.h"
#include "Eigen/LU"


using namespace Eigen;

//put this struct into global.h
struct SVDResult
{
    Eigen::Matrix3f U;
    Eigen::Matrix3f SIGMA;
    Eigen::Matrix3f V;
};

SVDResult SingularValueDecomposition3D(Matrix3f F)
{
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeThinU | Eigen::ComputeThinV);
    SVDResult result;
    Matrix3f tempU = svd.matrixU();
    Matrix3f tempV = svd.matrixV();

    Matrix3f tempSigma = svd.singularValues();

    //sorting
    if(tempU.determinant()<0)
    {
        tempU(2,0) *= -1;
        tempU(2,1) *= -1;
        tempU(2,2) *= -1;
        tempSigma(2,2) *= -1;
    }

    if(tempV.determinant()<0)
    {
        tempV(2,0) *= -1;
        tempV(2,1) *= -1;
        tempV(2,2) *= -1;
        tempSigma(2,2) *= -1;
    }

    if(tempSigma(0,0)<tempSigma(1,1))
    {
        float tempRecord = tempSigma(0,0);
        tempSigma(0,0) = tempSigma(1,1);
        tempSigma(1,1) = tempRecord;
    }

    result.U = tempU;
    result.V = tempV;
    result.SIGMA = tempSigma;

    return result;
}