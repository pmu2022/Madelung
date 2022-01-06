
#include <gtest/gtest.h>

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>

#include "spherical_harmonics.hpp"

/**
 * Test the CLM prefactor creation routines
 */
TEST(SphericalHarmonicsTest, Clm1) {

  int lmax = 3;

  auto *clm_local = new double[(lmax + 1) * (lmax + 2) / 2];

  calc_clm(&lmax, clm_local);

  // l = 0  m = 0
  const auto clm_ref1 = sqrt(1.0 / (4.0 * M_PI));
  EXPECT_NEAR(clm_ref1, clm_local[0], 1.0e-12);

  // l = 1  m = 0
  const auto clm_ref2 = sqrt(3.0 / (4.0 * M_PI));
  EXPECT_NEAR(clm_ref2, clm_local[1], 1.0e-12);

  delete[] clm_local;
}

/**
 * Test the CLM prefactor creation routines
 */
TEST(SphericalHarmonicsTest, Clm2) {

  /*
   *   0.28209479177387831
  0.48860251190292020
 -0.34549414947133572
  0.63078313050504031
 -0.25751613468212653
  0.12875806734106326
   */

  int lmax = 6;
  auto *clm_local = new double[(lmax + 1) * (lmax + 2) / 2];

  calc_clm(&lmax, clm_local);

  EXPECT_NEAR(0.28209479177387831, clm_local[0], 1.0e-12);
  EXPECT_NEAR(0.48860251190292020, clm_local[1], 1.0e-12);
  EXPECT_NEAR(-0.34549414947133572, clm_local[2], 1.0e-12);
  EXPECT_NEAR(0.63078313050504031, clm_local[3], 1.0e-12);
  EXPECT_NEAR(-0.25751613468212653, clm_local[4], 1.0e-12);
  EXPECT_NEAR(0.12875806734106326, clm_local[5], 1.0e-12);

  delete[] clm_local;
}

/**
 * Test the spherical harmonics `sph_harm_0` and `sph_harm_1`
 */
TEST(SphericalHarmonics, SphHarmTest1) {

  int lmax = 3;

  auto *ylm = new std::complex<double>[(lmax + 1) * (lmax + 1)];

  double x = 1.0;
  double y = 1.0;
  double z = 1.0;

  sph_harm_0(&x, &y, &z, &lmax, ylm);

  // l = 0  m = 0
  const auto real_ref00 = sqrt(1.0 / (4.0 * M_PI));
  const auto imag_ref00 = 0.0;

  EXPECT_NEAR(real_ref00, std::real(ylm[0]), 1.0e-12);
  EXPECT_NEAR(imag_ref00, std::imag(ylm[0]), 1.0e-12);

  std::vector<double> vec = {1.0, 1.0, 1.0};

  sph_harm_1(vec.data(), &lmax, ylm);

  // l = 0  m = 0
  EXPECT_NEAR(real_ref00, std::real(ylm[0]), 1.0e-12);
  EXPECT_NEAR(imag_ref00, std::imag(ylm[0]), 1.0e-12);


  delete[] ylm;

}

/**
 * Test spherical harmonics that are used for the calculation of the Madelung constant
 * and reduced Green's function
 */
TEST(SphericalHarmonics, SphHarmTest2) {

  int lmax = 3;

  auto *ylm = new std::complex<double>[(lmax + 1) * (lmax + 1)];

  std::vector<double> vec(3);

  vec[0] = -4.8830810748721376;
  vec[1] = -4.8830810748721376;
  vec[2] = -4.8830810748721376;

  sph_harm_1(vec.data(), &lmax, ylm);

  // (0.28209479177387831,0.0000000000000000)
  EXPECT_NEAR(0.28209479177387831, std::real(ylm[0]), 1.0e-10);
  EXPECT_NEAR(0.0, std::imag(ylm[0]), 1.0e-10);

  // (-0.19947114020071646,0.19947114020071652)
  EXPECT_NEAR(-0.19947114020071646, std::real(ylm[1]), 1.0e-10);
  EXPECT_NEAR(0.19947114020071652, std::imag(ylm[1]), 1.0e-10);

  // (-0.28209479177387831,-0.0000000000000000)
  EXPECT_NEAR(-0.28209479177387831, std::real(ylm[2]), 1.0e-10);
  EXPECT_NEAR(0.0, std::imag(ylm[2]), 1.0e-10);

  // (0.19947114020071646,0.19947114020071652)
  EXPECT_NEAR(0.19947114020071646, std::real(ylm[3]), 1.0e-10);
  EXPECT_NEAR(0.19947114020071652, std::imag(ylm[3]), 1.0e-10);

  delete[] ylm;

}

