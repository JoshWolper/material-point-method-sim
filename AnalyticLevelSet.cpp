//
// Created by fanfu on 11/24/17.
//

#include "AnalyticLevelSet.h"

bool AnalyticLevelSet::inside(const Vector3f& X) const {
    return signedDistance(X) <= 0;
}

Sphere::Sphere(const Vector3f &center, const float &radius) : center(center), radius(radius){}

float Sphere::signedDistance(const Vector3f& X) const {
    return (X-center).norm() - radius;
}

bool Sphere::inside(const Vector3f &X) const {
    return (X-center).squaredNorm() <= radius*radius;
}

Vector3f Sphere::normal(const Vector3f &X) const {
    Vector3f outward_normal = X - center;
    if (outward_normal.squaredNorm() < (float)1e-7)
        return Vector3f::Zero();
    return (outward_normal).normalized();
}