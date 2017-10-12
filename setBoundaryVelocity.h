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
    for (int i = 0; i<thickness; i++){
        for (int j = 0; j < Y; j++){
            for (int k = 0; k < Z; k++){
                index = i*Y*Z + j*Z + k;
                GridAttr thisGrid = gridAttrs[index];
                if (thisGrid.velG[0] < 0){

                }
            }
        }
    }
}

#endif //MPM_SETBOUNDARYVELOCITY_H
