/**
 *
 * Regression tests to ensure that the results are the same all the time
 *

 */


#include <gtest/gtest.h>

#define _USE_MATH_DEFINES

#include <cmath>
#include <cassert>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <tuple>

#include "multipole_madelung.h"

#include "MultipoleMadelung.hpp"

/**
 * Test for simple structure with two atoms
 */
TEST(RegressionTestSuite, Structure1) {

  int num_local_atoms = 2;
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

  EXPECT_NEAR(-0.2837297479481e+01, madelung.getMadSum(0, 0), 1e-12);
  EXPECT_NEAR(-0.8019359700280e+00, madelung.getMadSum(1, 0), 1e-12);
  EXPECT_NEAR(-0.8019359700280e+00, madelung.getMadSum(0, 1), 1e-12);
  EXPECT_NEAR(-0.2837297479481e+01, madelung.getMadSum(1, 1), 1e-12);

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

