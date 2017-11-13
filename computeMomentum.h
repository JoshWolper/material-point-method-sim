//
// Created by ziyinqu on 11/12/17.
//

#ifndef MPM_COMPUTEMOMENTUM_H
#define MPM_COMPUTEMOMENTUM_H

#include "global.h"
#include "Eigen/Eigen"

Vector3f computeParticleMomentum(const std::vector<Particle> particles){
    Vector3f Lp = Vector3f::Zero();
    for (int i = 0; i < particles.size(); i++){
        Lp += particles[i].massP * particles[i].velP;
    }
    return Lp;
}

Vector3f computeGridMomentum(const std::vector<GridAttr> gridAttrs, bool flag){
    Vector3f Lg = Vector3f::Zero();
    if (flag == false){
        for (int i = 0; i < gridAttrs.size(); i++){
            Lg += gridAttrs[i].massG * gridAttrs[i].velGn;
        }
    }
    else{
        for (int i = 0; i < gridAttrs.size(); i++){
            Lg += gridAttrs[i].massG * gridAttrs[i].velG;
        }
    }
    return Lg;
}
#endif //MPM_COMPUTEMOMENTUM_H
