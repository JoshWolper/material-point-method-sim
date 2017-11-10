//
// Created by ziyinqu on 11/9/17.
//

#include <iostream>
#include "kernelTest.h"
#include "../global.h"
#include "../Eigen/Eigen"
#include "../interpolation.h"

using namespace Eigen;

void quadraticTest(){
    std::cout << "INFO: >>>>>>>>>>>>>>>>>>> Kernel Test Start <<<<<<<<<<<<<<<<<<<" << std::endl;
    Vector3f xp = Vector3f(0.321932, 0.452119, 0.434341);
    float h = 0.04;
    Vector3f index_space = xp / h;
    Vector3i base_node;
    Matrix3f wp;
    Matrix3f dwp;
    QuadraticInterpolation(index_space, base_node, wp, dwp);
    std::cout << "INFO: " << xp.transpose() << std::endl;
    std::cout << "INFO: " << base_node.transpose() << std::endl;
    for (int i = 0; i < 3; i++){
        float wi = wp(0,i);
        float dwi = dwp(0,i);
        int node_i = base_node(0) + i;
        for (int j = 0; j < 3; j++){
            float wj = wp(1,j);
            float dwj = dwp(1,j);
            int node_j = base_node(1) + j;
            for (int k = 0; k < 3; k++){
                float wk = wp(2,k);
                float dwk = dwp(2,k);
                int node_k = base_node(2) + k;
                Vector3i node = Vector3i(node_i, node_j, node_k);
                float weight = wi * wj * wk;
                Vector3f grad_weight = Vector3f(dwi*wj*wk/h, wi*dwj*wk/h, wi*wj*dwk/h);
                std::cout << "INFO: " << "ijk: " << node.transpose() << std::endl;
                std::cout << "INFO: " << "weight: " << weight << std::endl;
                std::cout << "INFO: " << "grad_weight: " << grad_weight.transpose() << std::endl;
            }
        }
    }

}
