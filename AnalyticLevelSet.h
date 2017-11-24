//
// Created by fanfu on 11/24/17.
//

#ifndef MPM_ANALYTICLEVELSET_H
#define MPM_ANALYTICLEVELSET_H

#include "Eigen/Eigen"
#include "global.h"

class AnalyticLevelSet {
public:
    virtual ~AnalyticLevelSet() = default;

    virtual float signedDistance(const Vector3f& X) const = 0;

    virtual bool inside(const Vector3f& X) const ;

    virtual Vector3f normal(const Vector3f& X) const = 0;
};

class Sphere: public AnalyticLevelSet{
public:
    float radius;
    Vector3f center;

    Sphere(const Vector3f& center, const float& radius);

    ~Sphere() override = default;

    bool inside(const Vector3f& X) const override;

    float signedDistance(const Vector3f& X) const override ;

    Vector3f normal(const Vector3f& X) const override;

};
#endif //MPM_ANALYTICLEVELSET_H
