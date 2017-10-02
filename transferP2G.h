//
// Created by ziyinqu on 10/1/17.
//

#ifndef MPM_TRANSFERP2G_H
#define MPM_TRANSFERP2G_H

#include "interpolation.h"
#include "global.h"

using namespace std;

void transferP2G(vector<Particle> &particles, vector<Grid> &grids);
void P2GKernel(Particle &par, Grid &grid, Matrix3f &wp, Matrix3f &dwp);

//void transferP2G(vector<Vector3f>& xp, vector<float>& mp, vector<Vector3f>& vp, vector<Matrix3f>& BMatrix, MatrixXd& gridM, MatrixXd& gridVn);
//void P2GKernel(Vector3f &xp, Vector3f &mp, Vector3f &vp,Matrix3f B);

#endif //MPM_TRANSFERP2G_H
