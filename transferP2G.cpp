//
// Created by ziyinqu on 10/1/17.
//

#include "transferP2G.h"


//void P2GKernel(Vector3f &xp, Vector3f &mp, Vector3f &vp,Matrix3f B)
//{
//    Vector3f newXp;
//    Vector3f
//}


//void transferP2G(vector<Vector3f>& xp, vector<float>& mp, vector<Vector3f>& vp, vector<Matrix3f>& BMatrix, MatrixXd& gridM, MatrixXd& gridVn){
//
//    int iterationNum = xp.size();
//    for(int i=0;i<iterationNum;i++)
//    {
//
//    }
//}

void transferP2G(vector<Particle> &particles,vector<Grid> &grids)
{
    int iterationNum = particles.size();
    for(int i =0;i<iterationNum;i++)
    {
        P2GKernel(particles[i],grids[i]);
    }
}

void P2GKernel(Particle &par,Grid &grid)
{
    Vector3i baseNode;
    Matrix3f inp;
    Matrix3f gradInp;

    QuadraticInterpolation(par.posP,baseNode,inp,gradInp);


}
