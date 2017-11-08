//
// Created by ziyinqu on 10/1/17.
//
#pragma once
#ifndef MPM_ADVECTION_H
#define MPM_ADVECTION_H

#include "global.h"

using namespace std;

void addGravity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, Vector3f gravity){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        gridAttrs[index].force = gridAttrs[index].force + gridAttrs[index].massG * gravity ;

    }
}

void updateGridvelocity(vector<GridAttr>& gridAttrs, vector<int> active_nodes, float dt){
    for (int i = 0; i < active_nodes.size(); i++){
        int index = active_nodes[i];
        Vector3f test = gridAttrs[index].force / gridAttrs[index].massG;
        gridAttrs[index].velG = gridAttrs[index].velGn + dt * gridAttrs[index].force / gridAttrs[index].massG;
        //gridAttrs[index].velG = Vector3f(0,-1,0);
    }
}


#endif //MPM_ADVECTION_H
