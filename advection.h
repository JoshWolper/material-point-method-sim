//
// Created by ziyinqu on 10/1/17.
//

#ifndef MPM_ADVECTION_H
#define MPM_ADVECTION_H

#include "global.h"
#include "Eigen/Eigen"

using namespace std;

inline void addGravity(vector<GridAttr>& gridAttrs, vector<int>& active_nodes, Vector3f gravity){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        gridAttrs[index].force = gridAttrs[index].force + gridAttrs[index].massG * gravity;
    }
}

inline void updateGridvelocity(vector<GridAttr> gridAttrs, vector<int>& active_nodes, float dt){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        gridAttrs[index].velG = gridAttrs[index].velGn + dt * gridAttrs[index].force / gridAttrs[index].massG;
    }
}

#endif //MPM_ADVECTION_H
