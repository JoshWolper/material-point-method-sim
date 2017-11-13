//
// Created by ziyinqu on 11/9/17.
//

#include "derivativeTest.h"
#include "../Eigen/Eigen"
#include "../constitutiveModel.h"
#include <iostream>

using namespace Eigen;

void derivativeTest()
{
    using T = double;
    using TV = Vector3d;
    using TM = Matrix3d;
    constexpr int dim = 3;

    TM F0 = TM::Identity();
    F0 *= 10;
    TM dF = TM::Random(dim,dim);
    std::cout << "F0 \n" << F0 << "\ndF = \n" << dF << std::endl;
    for(int i = 2; i < 7; ++i) {
        T epsilon = (T)std::pow((T)10, -(T)i);
        TM F1 = F0 + epsilon * dF;
        TM piola0, piola1;
        T energy0, energy1;
        //corotatedPiolaDouble(F0, energy0, piola0);
        //corotatedPiolaDouble(F1, energy1, piola1);
        neoHookeanPiolaDouble(F0, energy0, piola0);
        neoHookeanPiolaDouble(F1, energy1, piola1);
        //stVernantPiolaDouble(F0, energy0, piola0);
        //stVernantPiolaDouble(F1, energy1, piola1);
        T contraction = 0;
        for(int row = 0; row < dim; ++row)
            for(int col = 0; col < dim; ++col)
                contraction += (piola0(row,col) + piola1(row,col)) * dF(row,col);
        contraction *= epsilon / 2;

        T dif = energy1 - energy0 - contraction;

        std::cout << "epsilon = " << epsilon << "\tenergy1 " << energy1 << "\tenergy0 = " << energy0 << "\tenergy1 - energy0 = " << energy1 - energy0 << "\tcontraction = " << contraction << "\tenergy1 - energy0 - contraction = " << dif << std::endl << std::flush;
    }
}