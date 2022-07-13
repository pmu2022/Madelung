/**
 *
 * Regression tests to ensure that the results are the same all the time
 *
 *             Madelung Matrix -- returns M_{0,0}(1:Na,i)
 *
 *             DL matrix   -- returns alat^l * DL(1:Na,1:jmax,i)
 *
 *             DL Factor   -- returns the prefactor of DL matrix.
 *                            The Madelung matrix is a product of
 *                            this factor and DL matrix.
 *
 */

#include <gtest/gtest.h>

#include <complex>
#include <iostream>
#include <vector>

#include "MultipoleMadelung.hpp"
#include "common.hpp"
#include "multipole_madelung.h"

/**
 * Tests for the development of the Madelung potential
 *
 * Supercell structure is tested
 *
 */
TEST(RegressionSupercellTestSuite, RegStructure3) {
  int num_local_atoms = 2;
  int num_atoms = 2;
  int gindex[] = {1, 2};
  int lmax_rho = 0;
  int lmax_pot = 0;
  int iprint = -10;

  double scale = 16;
  double a = 5.3821038770852745;

  double bravais[] = {
      1.0 * a * scale, 0.0, 0.0,
      0.0, 1.0 * a * scale, 0.0,
      0.0, 0.0, 1.0 * a * scale,
  };

  double position[] = {0.0, 0.0, 0.0, 0.5 * a, 0.5 * a, 0.5 * a};

  initMadelung(num_local_atoms, num_atoms, gindex, lmax_rho, lmax_pot, bravais,
               position, 10);

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

