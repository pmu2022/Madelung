/**
 *
 * Regression tests to ensure that the results are the same all the time
 *
 */


#include <gtest/gtest.h>

#include <complex>
#include <iostream>
#include <vector>

#include "multipole_madelung.h"

#include "MultipoleMadelung.hpp"

/**
 * Test for simple structure with two atoms
 */
TEST(RegressionTestSuite, Structure1) {

  int num_atoms = 2;
  std::vector<int> gindex{0, 1};
  int lmax = 3;

  Matrix<double> lattice(3, 3);

  lattice = 0;
  lattice(0, 0) = 1.0;
  lattice(1, 1) = 1.0;
  lattice(2, 2) = 1.0;

  Matrix<double> position(3, num_atoms);

  position(0, 0) = 0.0;
  position(1, 0) = 0.0;
  position(2, 0) = 0.0;

  position(0, 1) = 0.5;
  position(1, 1) = 0.5;
  position(2, 1) = 0.5;

  lsms::MultipoleMadelung madelung(lattice, position, lmax, gindex);

  // Test of all madelung values
  EXPECT_NEAR(-0.2837297479481e+01, madelung.getMadSum(0, 0), 1e-12);
  EXPECT_NEAR(-0.8019359700280e+00, madelung.getMadSum(1, 0), 1e-12);
  EXPECT_NEAR(-0.8019359700280e+00, madelung.getMadSum(0, 1), 1e-12);
  EXPECT_NEAR(-0.2837297479481e+01, madelung.getMadSum(1, 1), 1e-12);

  // Test of all dl matrix values
  EXPECT_NEAR(-10.057957687339862, std::real(madelung.getDlMatrix(0, 0, 0)), 1e-12);
  EXPECT_NEAR(0.0, std::imag(madelung.getDlMatrix(0, 0, 0)), 1e-12);

  // Test of all dl factors
  EXPECT_NEAR( 0.28209479177387842, madelung.getDlFactor(0, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959593E-002, madelung.getDlFactor(0, 1), 1e-12);
  EXPECT_NEAR( -9.4031597257959454E-002, madelung.getDlFactor(0, 2), 1e-12);
  EXPECT_NEAR(1.8806319451591908E-002, madelung.getDlFactor(0, 3), 1e-12);
  EXPECT_NEAR(-1.8806319451591908E-002, madelung.getDlFactor(0, 4), 1e-12);

}


/**
 * Tests for the development of the Madelung potential
 */
TEST(RegressionTestSuite, RegStructure1) {

  int num_local_atoms = 2;
  int num_atoms = 2;
  int gindex[] = {1, 2};
  int lmax_rho = 3;
  int lmax_pot = 3;
  int iprint = 10;

  double bravais[] = {
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
  };

  double position[] = {
      0.0, 0.0, 0.0,
      0.5, 0.5, 0.5
  };

  initMadelung(num_local_atoms,
               num_atoms,
               gindex,
               lmax_rho,
               lmax_pot,
               bravais,
               position,
               10);

  /*
   * Madelung Matrix(j,i) for atom     1    1: -0.2837297479481D+01
   * Madelung Matrix(j,i) for atom     2    1: -0.8019359700280D+00
   * Sum over j of Madelung Matrix(j,i) for atom     1:  -0.3639233449509D+01
   * Madelung Matrix(j,i) for atom     1    2: -0.8019359700280D+00
   * Madelung Matrix(j,i) for atom     2    2: -0.2837297479481D+01
   * Sum over j of Madelung Matrix(j,i) for atom     2:  -0.3639233449509D+01
   */

  printMadelungMatrix(&iprint);
  endMadelung();

}

