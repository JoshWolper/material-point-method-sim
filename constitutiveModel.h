//
// Created by ziyinqu on 11/9/17.
//
#pragma once
#ifndef MPM_CONSTITUTIVEMODEL_H
#define MPM_CONSTITUTIVEMODEL_H

#include "interpolation.h"
#include "SVD.h"

void corotatedPiola(Matrix3f defGrad, Matrix3f& piola);

void corotatedPiolaDouble(Matrix3d defGrad, double& energy, Matrix3d& piola);

void neoHookeanPiola(Matrix3f defGrad, Matrix3f& piola);

void neoHookeanPiolaDouble(Matrix3d defGrad, double& energy, Matrix3d& piola);

void stVernantPiola(Matrix3f defGrad, Matrix3f& piola);

void stVernantPiolaDouble(Matrix3d defGrad, double& energy, Matrix3d& piola);

void snowPiola(Matrix3f defGrad, Matrix3f Fp, Matrix3f Fe, Matrix3f& piola);
#endif //MPM_CONSTITUTIVEMODEL_H
