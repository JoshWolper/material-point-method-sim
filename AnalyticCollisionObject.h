//
// Created by fanfu on 11/24/17.
//

#ifndef MPM_ANALYTICCOLLISIONOBJECT_H
#define MPM_ANALYTICCOLLISIONOBJECT_H

#include "AnalyticLevelSet.h"
#include "global.h"
#include "Eigen/Eigen"

void SphereGridCollision(std::vector<GridAttr>& gridAttr, const std::vector<int>& active_nodes, const Vector3f& center, float radius, float friction, float dx);

void SphereParticleCollision(std::vector<Particle>& particles, const Vector3f& center, float friction, float radius);
#endif //MPM_ANALYTICCOLLISIONOBJECT_H
