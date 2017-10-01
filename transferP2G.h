//
// Created by ziyinqu on 10/1/17.
//

#ifndef MPM_TRANSFERP2G_H
#define MPM_TRANSFERP2G_H

#include "Eigen/Eigen"

using namespace Eigen;
using namespace std;

void transferP2G(vector<Vector3f>& xp, vector<float>& mp, vector<Vector3f>& vp, MatrixXd& gridM, MatrixXd& gridVn);

#endif //MPM_TRANSFERP2G_H
