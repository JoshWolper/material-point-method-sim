//
// Created by fanfu on 11/24/17.
//

#include "AnalyticCollisionObject.h"

#define STICKY true

void SphereGridCollision(std::vector<GridAttr>& gridAttr, const std::vector<int>& active_nodes, const Vector3f& center, float radius, float friction, float dx){
    // initialize sphere level set
    Sphere SphereLevelSet(center, radius);

    Vector3f xi = Vector3f::Zero();
    Vector3f n = Vector3f::Zero();

    for (int i = 0; i < active_nodes.size(); i++){

        int index = active_nodes[i];
        xi = dx * gridAttr[index].Xi;

        if (SphereLevelSet.inside(xi)){
            if (STICKY){
                gridAttr[index].velG = Vector3f::Zero();
            }
            else{
                n = SphereLevelSet.normal(xi);
                float dot = n.dot(gridAttr[index].velG);
                if (dot < 0){
                    gridAttr[index].velG -=  dot * n;
                    if (friction != 0){
                        if (-dot * friction < gridAttr[index].velG.norm())
                            gridAttr[index].velG += dot * friction * gridAttr[index].velG.normalized();
                        else
                            gridAttr[index].velG = Vector3f::Zero();
                    }
                }
            }


        }
    }
}

void SphereParticleCollision(std::vector<Particle>& particles, const Vector3f& center, float friction, float radius){
    // initialize sphere level set
    Sphere SphereLevelSet(center, radius);

    Vector3f xp = Vector3f::Zero();
    Vector3f n = Vector3f::Zero();

    for (int i = 0; i < particles.size(); i++){
        if (SphereLevelSet.inside(xp)){
            if (STICKY){
                particles[i].velP = Vector3f::Zero();
            }
            else{
                n = SphereLevelSet.normal(xp);
                float dot = n.dot(particles[i].velP);
                if (dot < 0){
                    particles[i].velP -=  dot * n;
                    if (friction != 0){
                        if (-dot * friction < particles[i].velP.norm())
                            particles[i].velP += dot * friction * particles[i].velP.normalized();
                        else
                            particles[i].velP = Vector3f::Zero();
                    }
                }
            }
        }
    }
}