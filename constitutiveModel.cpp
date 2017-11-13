//
// Created by ziyinqu on 11/9/17.
//

#include "constitutiveModel.h"
#include <cmath>

void corotatedPiolaDouble(Matrix3d defGrad, double& energy, Eigen::Matrix3d& piola){

    double E = 100;
    double nu = 0.3;

    double lambda = E * nu / (((double)1 + nu) * ((double)1 - (double)2 * nu));
    double mu = E / ((double)2 * ((double)1 + nu));

    //cout << "DefGrad: " << defGrad << endl;

    SVDResultDouble svdResult = SingularValueDecomposition3DDouble(defGrad);

    Matrix3d U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    //cout << "U: " << U << endl;
    //cout << "Sigma: " << sigma << endl;
    //cout << "V: " << V << endl;

    Matrix3d R = U * V.transpose();

    double J = defGrad.determinant();

    piola = (2 * mu * (defGrad - R)) + (lambda * (J-1) * J * (defGrad.inverse().transpose()));

    energy = mu * (defGrad - R).squaredNorm() + (lambda / 2) * (J - 1) * (J - 1);


    return;
}

void corotatedPiola(Matrix3f defGrad, Eigen::Matrix3f& piola){

    float E = 10;
    float nu = 0.3;

    float lambda = E * nu / (((float)1 + nu) * ((float)1 - (float)2 * nu));
    float mu = E / ((float)2 * ((float)1 + nu));

    //cout << "DefGrad: " << defGrad << endl;

    SVDResult svdResult = SingularValueDecomposition3D(defGrad);

    Matrix3f U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    //cout << "U: " << U << endl;
    //cout << "Sigma: " << sigma << endl;
    //cout << "V: " << V << endl;

    Matrix3f R = U * V.transpose();

    float J = defGrad.determinant();
    Matrix3f JFinvT;
//    JFinvT(0, 0) = defGrad(1, 1) * defGrad(2, 2) - defGrad(1, 2) * defGrad(2, 1);
//    JFinvT(0, 1) = defGrad(1, 2) * defGrad(2, 0) - defGrad(1, 0) * defGrad(2, 2);
//    JFinvT(0, 2) = defGrad(1, 0) * defGrad(2, 1) - defGrad(1, 1) * defGrad(2, 0);
//    JFinvT(1, 0) = defGrad(0, 2) * defGrad(2, 1) - defGrad(0, 1) * defGrad(2, 2);
//    JFinvT(1, 1) = defGrad(0, 0) * defGrad(2, 2) - defGrad(0, 2) * defGrad(2, 0);
//    JFinvT(1, 2) = defGrad(0, 1) * defGrad(2, 0) - defGrad(0, 0) * defGrad(2, 1);
//    JFinvT(2, 0) = defGrad(0, 1) * defGrad(1, 2) - defGrad(0, 2) * defGrad(1, 1);
//    JFinvT(2, 1) = defGrad(0, 2) * defGrad(1, 0) - defGrad(0, 0) * defGrad(1, 2);
//    JFinvT(2, 2) = defGrad(0, 0) * defGrad(1, 1) - defGrad(0, 1) * defGrad(1, 0);
    piola = (2 * mu * (defGrad - R)) + (lambda * (J-1) * J * (defGrad.transpose().inverse()));
//    piola = (2 * mu * (defGrad - R)) + (lambda * (J-1) * JFinvT);
    return;
}

void neoHookeanPiola(Matrix3f defGrad, Eigen::Matrix3f& piola){

    float E = 10;
    float nu = 0.3;

    float mu = E / ((float)2 * ((float)1 + nu));
    float lambda = E * nu / (((float)1 + nu) * ((float)1 - ((float)2*nu)));

    SVDResult svdResult = SingularValueDecomposition3D(defGrad);

    Matrix3f U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    Matrix3f R = U * V.transpose();

    float J = defGrad.determinant();

    piola = (mu * (defGrad - defGrad.transpose().inverse())) + (lambda * log(J) * defGrad.transpose().inverse());

    return;
}

void neoHookeanPiolaDouble(Matrix3d defGrad, double& energy, Eigen::Matrix3d& piola){

    double E = 500000;
    double nu = 0.3;

    double mu = E / ((double)2 * ((double)1 + nu));
    double lambda = E * nu / (((double)1 + nu) * ((double)1 - ((double)2*nu)));

    SVDResultDouble svdResult = SingularValueDecomposition3DDouble(defGrad);

    Matrix3d U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    Matrix3d R = U * V.transpose();

    double J = defGrad.determinant();

    piola = (mu * (defGrad - defGrad.transpose().inverse())) + (lambda * log(J) * defGrad.transpose().inverse());

    Matrix3d fTf = defGrad.transpose() * defGrad;
    energy = ((mu / (double)2) * (fTf.trace() - (double)3)) - (mu * log(J)) + ((lambda/(double)2) * log(J) * log(J));

    return;
}

void stVernantPiola(Matrix3f defGrad, Eigen::Matrix3f& piola){

    float E = 10;
    float nu = 0.3;

    float mu = E / ((float)2 * ((float)1 + nu));
    float lambda = E * nu / (((float)1 + nu) * ((float)1 - ((float)2*nu)));

    SVDResult svdResult = SingularValueDecomposition3D(defGrad);

    Matrix3f U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    Matrix3f R = U * V.transpose();

    float J = defGrad.determinant();

    Matrix3f logSigma = Matrix3f::Zero();
    logSigma(0,0) = log(sigma(0,0));
    logSigma(1,1) = log(sigma(1,1));
    logSigma(2,2) = log(sigma(2,2));

    Matrix3f piolaSingular = (2 * mu * logSigma * sigma.inverse()) + (lambda * logSigma.trace() * sigma.inverse()); //calculate singular value view of piola

    //now calculate actual P!
    piola = U * piolaSingular * V.transpose();

    return;
}

void stVernantPiolaDouble(Matrix3d defGrad, double& energy, Eigen::Matrix3d& piola){

    double E = 500000;
    double nu = 0.3;

    double mu = E / ((double)2 * ((double)1 + nu));
    double lambda = E * nu / (((double)1 + nu) * ((double)1 - ((double)2*nu)));

    SVDResultDouble svdResult = SingularValueDecomposition3DDouble(defGrad);

    Matrix3d U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    Matrix3d R = U * V.transpose();

    float J = defGrad.determinant();

    Matrix3d logSigma = Matrix3d::Zero();
    logSigma(0,0) = log(sigma(0,0));
    logSigma(1,1) = log(sigma(1,1));
    logSigma(2,2) = log(sigma(2,2));

    Matrix3d piolaSingular = (2 * mu * logSigma * sigma.inverse()) + (lambda * logSigma.trace() * sigma.inverse()); //calculate singular value view of piola

    //now calculate actual P!
    piola = U * piolaSingular * V.transpose();

    energy = (mu * (logSigma * logSigma).trace()) + ((lambda/(double)2) * (logSigma.trace()) * logSigma.trace());

    return;
}