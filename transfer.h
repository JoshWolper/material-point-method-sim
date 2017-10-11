//
// Created by ziyinqu on 10/2/17.
//

#ifndef MPM_TRANSFER_H
#define MPM_TRANSFER_H

#include "global.h"
#include "Eigen/Eigen"
#include "interpolation.h"

using namespace std;

void transferG2P(vector<Particle>& particles, vector<GridAttr>& gridAttrs, const GridInfo gridInfo, float dt, float alpha){
    int iterationNum = particles.size();
    Matrix3f wp;
    Matrix3f dwp;
    Vector3i baseNode;
    Vector3f vpic;
    Vector3f vflip;
    int H = gridInfo.H;
    int L = gridInfo.L;
    for(int iter =0; iter < iterationNum; iter++)
    {
        QuadraticInterpolation(particles[iter].posP, baseNode, wp, dwp);
        for (int i = 0; i < 3; i++){
            float wi = wp(0,i);
            int node_i = baseNode(0) + (i-1);
            for (int j = 0; j < 3; j++){
                float wij = wi*wp(1,j);
                int node_j = baseNode(1) + (j-1);
                for (int k = 0; k < 3; k++){
                    float wijk = wij*wp(2,k);
                    int node_k = baseNode(2) + (k-1);
                    int index = node_i*H*L + node_j*L + node_k;
                    // calculate PIC and FLIP part velocity
                    vpic += wijk * gridAttrs[index].velG;
                    vflip += wijk * (gridAttrs[index].velG - gridAttrs[index].velGn);
                }
            }
        }
        // update particles velocity and position
        particles[iter].velP = (1-alpha)*vpic + alpha*vflip;
        particles[iter].posP += dt * particles[iter].velP;
    }
}

void transferP2G(vector<Particle>& particles, vector<GridAttr>& gridAttrs, const GridInfo gridInfo, std::vector<int>& active_nodes)
{
    int iterationNum = particles.size();
    Matrix3f wp;
    Matrix3f dwp;
    Vector3i baseNode;
    int H = gridInfo.H;
    int L = gridInfo.L;
    for(int iter =0; iter < iterationNum; iter++)
    {
        QuadraticInterpolation(particles[iter].posP, baseNode, wp, dwp);
        for (int i = 0; i < 3; i++){
            float wi = wp(0,i);
            int node_i = baseNode(0) + (i-1);
            for (int j = 0; j < 3; j++){
                float wij = wi*wp(1,j);
                int node_j = baseNode(1) + (j-1);
                for (int k = 0; k < 3; k++){
                    float wijk = wij*wp(2,k);
                    int node_k = baseNode(2) + (k-1);
                    int index = node_i*H*L + node_j*L + node_k;
                    // grid mass transfer
                    gridAttrs[index].massG += particles[i].massP * wijk;
                    // grid velocity transfer
                    //TODO APIC transfer should apply here
                    gridAttrs[index].velGn += particles[i].massP * particles[i].velP;
                }
            }
        }
    }

    for (int iter = 0; iter < gridAttrs.size(); iter++){
        if (gridAttrs[iter].massG != 0){
            active_nodes.push_back(iter);
            gridAttrs[iter].velGn = gridAttrs[iter].velGn / gridAttrs[iter].massG ;
        }
    }
}
#endif //MPM_TRANSFER_H
