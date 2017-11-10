//
// Created by ziyinqu on 11/9/17.
//

#include "SVD.h"
#include "Eigen/Dense"
#include "Eigen/src/SVD/JacobiSVD.h"
#include "Eigen/LU"

SVDResult SingularValueDecomposition3D(Matrix3f F)
{
    Eigen::JacobiSVD<Eigen::Matrix3f> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SVDResult result;
    Matrix3f tempU = svd.matrixU();
    Matrix3f tempV = svd.matrixV();

    Vector3f singVals = svd.singularValues();
    Matrix3f tempSigma = Matrix3f::Zero();
    tempSigma(0,0) = singVals(0);
    tempSigma(1,1) = singVals(1);
    tempSigma(2,2) = singVals(2);

    //sorting
    if(tempU.determinant()<0)
    {
        tempU(0,2) *= -1;
        tempU(1,2) *= -1;
        tempU(2,2) *= -1;
        tempSigma(2,2) *= -1;
    }

    if(tempV.determinant()<0)
    {
        tempV(0,2) *= -1;
        tempV(1,2) *= -1;
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

SVDResultDouble SingularValueDecomposition3DDouble(Matrix3d F)
{
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    SVDResultDouble result;
    Matrix3d tempU = svd.matrixU();
    Matrix3d tempV = svd.matrixV();

    Vector3d singVals = svd.singularValues();
    Matrix3d tempSigma = Matrix3d::Zero();
    tempSigma(0,0) = singVals(0);
    tempSigma(1,1) = singVals(1);
    tempSigma(2,2) = singVals(2);

    //sorting
    if(tempU.determinant()<0)
    {
        tempU(0,2) *= -1;
        tempU(1,2) *= -1;
        tempU(2,2) *= -1;
        tempSigma(2,2) *= -1;
    }

    if(tempV.determinant()<0)
    {
        tempV(0,2) *= -1;
        tempV(1,2) *= -1;
        tempV(2,2) *= -1;
        tempSigma(2,2) *= -1;
    }

    if(tempSigma(0,0)<tempSigma(1,1))
    {
        double tempRecord = tempSigma(0,0);
        tempSigma(0,0) = tempSigma(1,1);
        tempSigma(1,1) = tempRecord;
    }

    result.U = tempU;
    result.V = tempV;
    result.SIGMA = tempSigma;
//    Matrix3d reconstruction = result.U * result.SIGMA * result.V.transpose();
//    Matrix3d dif = F - reconstruction;
//    std::cout << "F - reconstruction = \n" << dif << std::endl << std::flush;

    return result;
}