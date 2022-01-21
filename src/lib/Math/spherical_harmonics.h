//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <complex>

extern "C"
{

void calc_clm(int *lmax, double *clm_local);

void sph_harm_0(double *x, double *y, double *z, int *lmax, std::complex<double> *ylm);

void sph_harm_1(double *vec, int *lmax, std::complex<double> *ylm);

}


#endif //SPHERICAL_HARMONICS_H
