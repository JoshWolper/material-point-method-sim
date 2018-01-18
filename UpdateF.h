//
// Created by VElysianP on 11/7/2017.
//

#ifndef MPM_UPDATEF_H
#define MPM_UPDATEF_H

#endif //MPM_UPDATEF_H

#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/src/SVD/JacobiSVD.h"
#include "Eigen/LU"
#include "global.h"
#include "interpolation.h"
#include "SVD.h"

using namespace std;
using namespace Eigen;

#define DIMENSIONX 3
#define DIMENSIONY 3
#define DIMENSIONZ 3


//this is the kernel function that only process one Fp update
//it is different from the MATLAB version of update all the Fps together
//target: update the Fp inside the particle
void UpdateF(float timeStep, const GridInfo gridInfo, vector<GridAttr> gridAttrs, vector<Particle>& particles, int energyDensityFunction)
{
    int particleNum = particles.size();

    for(int loop = 0; loop < particleNum; loop++)
    {
        Vector3f particlePos = particles[loop].posP/gridInfo.dx;
        Vector3i baseNode;
        Matrix3f wp;
        Matrix3f dwp;
        QuadraticInterpolation(particlePos, baseNode, wp, dwp);

        Matrix3f defGrad = particles[loop].F;

        Matrix3f grad_vp = Matrix3f::Zero();

        for(int k = baseNode(2); k <= baseNode(2)+2 ; k++)
        {
            for(int j = baseNode(1); j <= baseNode(1)+2; j++)
            {
                for(int i = baseNode(0); i <= baseNode(0)+2;i++)
                {
                    int tempNodeIndex = i * gridInfo.H * gridInfo.L + j * gridInfo.L + k ;
                    Vector3i tempNodeLength = Vector3i(i,j,k) - baseNode;

                    Vector3f dwip = Vector3f((1/gridInfo.dx) * dwp(0,tempNodeLength(0))*wp(1,tempNodeLength(1))*wp(2,tempNodeLength(2)),
                                             (1/gridInfo.dx) * wp(0,tempNodeLength(0))*dwp(1,tempNodeLength(1))*wp(2,tempNodeLength(2)),
                                             (1/gridInfo.dx) * wp(0,tempNodeLength(0))*wp(1,tempNodeLength(1))*dwp(2,tempNodeLength(2)));

                    Vector3f tempVel = gridAttrs[tempNodeIndex].velG;
                    //TODO: bottleneck
                    Matrix3f test;

                    test(0,0) = tempVel(0) * dwip(0);
                    test(0,1) = tempVel(0) * dwip(1);
                    test(0,2) = tempVel(0) * dwip(2);
                    test(1,0) = tempVel(1) * dwip(0);
                    test(1,1) = tempVel(1) * dwip(1);
                    test(1,2) = tempVel(1) * dwip(2);
                    test(2,0) = tempVel(2) * dwip(0);
                    test(2,1) = tempVel(2) * dwip(1);
                    test(2,2) = tempVel(2) * dwip(2);
                    //grad_vp += tempVel * dwip.transpose();
                    grad_vp += test;
                }
            }
        }

        Matrix3f newdefGrad;
        Matrix3f newFe, newFp, Ue, sigmae, Ve;
        newdefGrad = defGrad + timeStep * grad_vp * defGrad;
        if (energyDensityFunction == 3){
            float thetaC = 2.5e-2;
            float thetaS = 5.5e-3;

            Matrix3f Fe = particles[loop].Fe;
            newFe = Fe + timeStep * grad_vp * Fe;
            SVDResult svdResult = SingularValueDecomposition3D(newFe);

            Matrix3f Ue, sigmae, Ve;
            Ue = svdResult.U;
            sigmae = svdResult.SIGMA;
            Ve= svdResult.V;

            for (int i = 0; i < 3; i++){
                sigmae(i,i) = max(1-thetaC, min(sigmae(i,i), 1+thetaS));
            }

            particles[loop].Fe = Ue * sigmae * Ve.transpose();
            particles[loop].Fp = Ve * sigmae.inverse() * Ue.transpose() * newdefGrad;

        }
        particles[loop].F = newdefGrad;
    }

}