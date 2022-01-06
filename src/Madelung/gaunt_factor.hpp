//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef MADELUNG_GAUNT_FACTOR_HPP
#define MADELUNG_GAUNT_FACTOR_HPP

#include "Matrix.hpp"
#include "Array3d.hpp"

extern "C"
{
void gaunt_factor(int *lmax, double *cgnt, int *kj3, int *nj3);
}

double get_gaunt_factor(Array3d<double> &cgnt,
                        Matrix<int> &nj3,
                        Array3d<int> &kj3,
                        int k1,
                        int k2,
                        int k3);



#endif //MADELUNG_GAUNT_FACTOR_HPP
