//
// Created by ziyinqu on 10/11/17.
//

#ifndef MPM_SETBOUNDARYVELOCITY_H
#define MPM_SETBOUNDARYVELOCITY_H

#include "global.h"

using namespace std;

void setBoundaryVelocity(vector<GridAttr>& gridAttrs, GridInfo gridInfo){
    int thickness = 2;
    int index = 0;
    int X = gridInfo.W;
    int Y = gridInfo.H;
    int Z = gridInfo.L;
    // boundary velocity for x-direction wall
    for (int i = 0; i<thickness; i++){
        for (int j = 0; j < Y; j++){
            for (int k = 0; k < Z; k++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[0] < 0){
                    gridAttrs[index].velG[0] = 0;
                }
            }
        }
    }
    for (int i = X-thickness; i<X; i++){
        for (int j = 0; j < Y; j++){
            for (int k = 0; k < Z; k++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[0] > 0){
                    gridAttrs[index].velG[0] = 0;
                }
            }
        }
    }
    // boudanry velocity for y-direction wall
    for (int j = 0; j<thickness; j++){
        for (int i = 0; i < X; i++){
            for (int k = 0; k < Z; k++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[1] < 0){
                    gridAttrs[index].velG[1] = 0;
                }
            }
        }
    }
    for (int j = Y-thickness; j<Y; j++){
        for (int i = 0; i < X; i++){
            for (int k = 0; k < Z; k++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[1] > 0){
                    gridAttrs[index].velG[1] = 0;
                }
            }
        }
    }
    // boudanry velocity for y-direction wall
    for (int k = 0; k<thickness; k++){
        for (int i = 0; i < X; i++){
            for (int j = 0; j < Y; j++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[2] < 0){
                    gridAttrs[index].velG[2] = 0;
                }
            }
        }
    }
    for (int k = Z-thickness; k<Z; k++){
        for (int i = 0; i < X; i++){
            for (int j = 0; j < Y; j++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[2] > 0){
                    gridAttrs[index].velG[2] = 0;
                }
            }
        }
    }
}

#endif //MPM_SETBOUNDARYVELOCITY_H
