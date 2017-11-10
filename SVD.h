//
// Created by VElysianP on 11/7/2017.
//
#pragma once
#ifndef MPM_SVD_H
#define MPM_SVD_H

#include "Eigen/Eigen"

#include "global.h"

using namespace Eigen;

SVDResult SingularValueDecomposition3D(Matrix3f F);

SVDResultDouble SingularValueDecomposition3DDouble(Matrix3d F);

#endif //MPM_SVD_H
