//
// Created by ziyinqu on 10/2/17.
//
#pragma once
#ifndef MPM_TRANSFER_H
#define MPM_TRANSFER_H

#include "global.h"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "interpolation.h"
#include "SVD.h"
#include <math.h>

#define USEAPIC false

using namespace std;

void transferG2P(vector<Particle>& particles, vector<GridAttr>& gridAttrs, const GridInfo gridInfo, float dt, float alpha){

    int iterationNum = particles.size();
    Matrix3f wp = Matrix3f::Zero();
    Matrix3f dwp = Matrix3f::Zero();
    Vector3i baseNode = Vector3i::Zero();
    int H = gridInfo.H;
    int L = gridInfo.L;

    for(int iter =0; iter < iterationNum; iter++) {
        Vector3f vpic = Vector3f::Zero();
        Vector3f vflip = Vector3f::Zero();
        Vector3f index_Space = particles[iter].posP/gridInfo.dx;
        QuadraticInterpolation(index_Space, baseNode, wp, dwp);
        for (int i = 0; i < 3; i++) {
            float wi = wp(0, i);
            int node_i = baseNode(0) + i;
            for (int j = 0; j < 3; j++) {
                float wij = wi * wp(1, j);
                int node_j = baseNode(1) + j;
                for (int k = 0; k < 3; k++) {
                    float wijk = wij * wp(2, k);
                    int node_k = baseNode(2) + k;
                    int index = node_i * H * L + node_j * L + node_k;

                    Vector3f gridIndex = Vector3f(node_i,node_j,node_k);
                    particles[iter].BP += wijk * gridAttrs[index].velG * (gridIndex.transpose() - index_Space.transpose());

                    // calculate PIC and FLIP part velocity
                    vpic += wijk * gridAttrs[index].velG;
                    vflip += wijk * (gridAttrs[index].velG - gridAttrs[index].velGn);
                }
            }
        }
        // update particles velocity and position

        if(!USEAPIC)
        {
            particles[iter].velP = (1 - alpha) * vpic + alpha * vflip;
        }
        else
        {
            particles[iter].velP = vpic;
        }

        particles[iter].posP += dt * particles[iter].velP;
    }
}

void transferP2G(vector<Particle>& particles, vector<GridAttr>& gridAttrs, const GridInfo gridInfo, std::vector<int>& active_nodes)
{
    int iterationNum = particles.size();
    Matrix3f wp = Matrix3f::Zero();
    Matrix3f dwp = Matrix3f::Zero();
    Vector3i baseNode = Vector3i::Zero();
    int H = gridInfo.H;
    int L = gridInfo.L;
    for(int iter =0; iter < iterationNum; iter++) {
        Vector3f index_Space = particles[iter].posP / gridInfo.dx;
        QuadraticInterpolation(index_Space, baseNode, wp, dwp);
        for (int i = 0; i < 3; i++) {
            float wi = wp(0, i);
            int node_i = baseNode(0) + i;
            for (int j = 0; j < 3; j++) {
                float wij = wi * wp(1, j);
                int node_j = baseNode(1) + j;
                for (int k = 0; k < 3; k++) {
                    float wijk = wij * wp(2, k);
                    int node_k = baseNode(2) + k;
                    int index = node_i * H * L + node_j * L + node_k;
                    // grid mass transfer
                    gridAttrs[index].massG += wijk * particles[iter].massP;

                    // calculate APIC things
                    // grid velocity transfer
                    //TODO APIC transfer should apply here
                    Vector3f gridNode = Vector3f(node_i, node_j, node_k);
                    Vector3f index_Space = particles[iter].posP / gridInfo.dx;
                    Vector3f plus = 4 * particles[i].BP * (gridNode - index_Space);

                    if (USEAPIC) {
                        gridAttrs[index].velGn += wijk * particles[i].massP * (particles[i].velP + plus);
                    } else {
                        gridAttrs[index].velGn += wijk * particles[i].massP * particles[i].velP;

                    }
                }
            }
        }
    }

    for (int iter = 0; iter < gridAttrs.size(); iter++){
        if (gridAttrs[iter].massG != 0){
            active_nodes.push_back(iter);
            gridAttrs[iter].velGn = gridAttrs[iter].velGn / gridAttrs[iter].massG ;
        }
        else{
            gridAttrs[iter].velGn = Vector3f::Zero();
        }
    }

}

void corotatedPiola(const Matrix3d& defGrad, double& energy, Eigen::Matrix3d& piola){

    double E = 10;
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

    piola = (2 * mu * (defGrad - R)) + (lambda * (J-1) * J * (defGrad.transpose().inverse()));



    return;
}

void neoHookeanPiola(Matrix3f defGrad, Eigen::Matrix3f& piola){

    float E = 500000;
    float nu = 0.3;

    float mu = E / (2 * (1 + nu));
    float lambda = E * nu / ((1 + nu) * (1 - (2*nu)));

    SVDResult svdResult = SingularValueDecomposition3D(defGrad);

    Matrix3f U, sigma, V;
    U = svdResult.U;
    sigma = svdResult.SIGMA;
    V = svdResult.V;

    Matrix3f R = U * V.transpose();

    float J = defGrad.determinant();

    piola = (mu * (defGrad - defGrad.transpose())) + (lambda * log(J) * defGrad.transpose().inverse());

    return;
}

void stVernantPiola(Matrix3f defGrad, Eigen::Matrix3f& piola){

    float E = 500000;
    float nu = 0.3;

    float mu = E / (2 * (1 + nu));
    float lambda = E * nu / ((1 + nu) * (1 - (2*nu)));

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

void addGridForces(vector<GridAttr>& gridAttrs, vector<Particle>& particles, GridInfo gridInfo, int energyDensityFunction) {
    //for (int i = 0; i < gridInfo.gridSize; i++){
    //    gridAttrs[i].force = Vector3f::Zero();
    //    cout << "Force: " << gridAttrs[i].force.transpose() << endl;
    //}

    //for each active grid node
    for (int i = 0; i < particles.size(); i++){
        //calculate force update
        float volume = particles[i].volumeP;
        Matrix3f defGrad = particles[i].F;
        float h = gridInfo.dx;

        //Calculate Piola Kirchoff Stress
        Matrix3f piola = Matrix3f::Ones();

        switch(energyDensityFunction){
            case 0: corotatedPiola(defGrad, piola);
            case 1: neoHookeanPiola(defGrad, piola);
            case 2: stVernantPiola(defGrad, piola);
            default: corotatedPiola(defGrad, piola); //default should just be the corotated model
        }

        Vector3f pos = particles[i].posP / gridInfo.dx;
        Vector3i baseNode;
        Matrix3f wp = Matrix3f::Ones();
        Matrix3f dwp = Matrix3f::Ones();

        QuadraticInterpolation(pos, baseNode, wp, dwp); //get Wp and gradWp

        //Loop through each of the 27 grid nodes
        for(int r = 0; r < 3; r++){
            for(int s = 0; s < 3; s++){
                for(int t = 0; t < 3; t++) {

                    Vector3f gradWip;
                    gradWip(0) = (1/h) * dwp(0,r) * wp(1, s) * wp(2, t); //calculate each component of gradWip
                    gradWip(1) = (1/h) * wp(0,r) * dwp(1, s) * wp(2, t);
                    gradWip(2) = (1/h) * wp(0,r) * wp(1, s) * dwp(2, t);

                    //cout << "Piola: " << piola << endl;
                    //cout << "defGrad: " << defGrad << endl;
                    //cout << "GradWip: " << gradWip << endl;

                    Vector3f f_i = -1 * volume * piola * defGrad.transpose() * gradWip; //calc force update

                    //cout << "f_i: " << gradWip.transpose() << endl;

                    int x = baseNode(0) + r; //calculate the indeces of the node we're acting on
                    int y = baseNode(1) + s;
                    int z = baseNode(2) + t;
                    int index = (x * gridInfo.H * gridInfo.L) + (y * gridInfo.L) + z;
                    gridAttrs[index].force = gridAttrs[index].force + f_i; //add the force

                }
            }
        }

    }
    /*
    for(int i = 0; i < gridInfo.gridSize; i++){
        if(std::abs(gridAttrs[i].force(0)) > 1e-10 || std::abs(gridAttrs[i].force(1)) > 1e-10 || std::abs(gridAttrs[i].force(2)) > 1e-10){
            cout << "Force: " << gridAttrs[i].force.transpose() << endl;
            //cout << flush;
        }
    }*/
}


#endif //MPM_TRANSFER_H
