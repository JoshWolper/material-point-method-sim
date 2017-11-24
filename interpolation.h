//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_INTERPOLATION_H
#define MPM_INTERPOLATION_H

#include "Eigen/Eigen"

Eigen::Vector3f calcWeights(float index_space, int& baseNode);

Eigen::Vector3f calcGradWeights(float index_space, int baseNode);

void QuadraticInterpolation(Eigen::Vector3f particlePos, Eigen::Vector3i& baseNode, Eigen::Matrix3f& wp, Eigen::Matrix3f& dwp);

void CubicInterpolation(Eigen::Vector3f particlePos, Eigen::Vector3i& baseNode, Eigen::Matrix4f& wp, Eigen::Matrix4f& dwp);
#endif //MPM_INTERPOLATION_H
