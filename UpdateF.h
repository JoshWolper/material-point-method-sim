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

using namespace std;
using namespace Eigen;

#define DIMENSIONX 3
#define DIMENSIONY 3
#define DIMENSIONZ 3


//this is the kernel function that only process one Fp update
//it is different from the MATLAB version of update all the Fps together
//target: update the Fp inside the particle
void UpdateF(int timeStep, const GridInfo gridInfo, vector<GridAttr> gridAttrs, vector<Particle>& particles)
{
    int particleNum = particles.size();

    for(int loop = 0; loop < particleNum; loop++)
    {
        Vector3f particlePos = particles[loop].posP/gridInfo.dx;
        Vector3i baseNode;
        Matrix3f wp;
        Matrix3f dwp;
        QuadraticInterpolation(particlePos, baseNode, wp, dwp);

        Matrix3f Fp = particles[loop].F;

        Matrix3f grad_vp = Matrix3f::Zero();

        for(int k = baseNode(2); k <= baseNode(2)+2 ; k++)
        {
            for(int j = baseNode(1); j <= baseNode(1)+2; j++)
            {
                for(int i = baseNode(0); i <= baseNode(0)+2;i++)
                {
                    int tempNodeIndex = i * gridInfo.H * gridInfo.L + j * gridInfo.L + k ;
                    Vector3i tempNodeLength = Vector3i(i,j,k) - baseNode;

                    Vector3f dwip = Vector3f(dwp(0,tempNodeLength(0))*wp(1,tempNodeLength(1))*wp(2,tempNodeLength(2)),
                                             wp(0,tempNodeLength(0))*dwp(1,tempNodeLength(1))*wp(2,tempNodeLength(2)),
                                             wp(0,tempNodeLength(0))*wp(1,tempNodeLength(1))*dwp(2,tempNodeLength(2)));

                    Vector3f tempVel = gridAttrs[tempNodeIndex].velG;
                    grad_vp += tempVel * dwip.transpose();
                }
            }
        }

        Matrix3f newFp = Fp + timeStep * grad_vp * Fp;
        particles[loop].F = newFp;
    }

}