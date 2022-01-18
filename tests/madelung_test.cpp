/**
 *
 * This tests cover all calculation that are nesseary for obtaining:
 *
 *         Madelung Matrix -- returns M_{0,0}(1:Na,i)
 *
 *             DL matrix   -- returns alat^l * DL(1:Na,1:jmax,i)
 *
 *             DL Factor   -- returns the prefactor of DL matrix.
 *                            The Madelung matrix is a product of
 *                            this factor and DL matrix.
 *
 */

#include <gtest/gtest.h>

#define _USE_MATH_DEFINES

#include <cassert>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

#include "MultipoleMadelung.hpp"
#include "lattice_utils.hpp"
#include "madelung.hpp"
#include "multipole_madelung.h"



/**
 * Test the scaling factors
 */
TEST(MadelungTests, ScalingFactor1) {
  int lmax = 3;

  Matrix<double> bravais(3, 3);
  bravais = 0.0;

  bravais(0, 0) = 1.0;
  bravais(1, 1) = 1.0;
  bravais(2, 2) = 1.0;

  auto scaling_factor = lsms::scaling_factor(bravais, lmax);
  EXPECT_NEAR(scaling_factor, 0.81915494309189596, 1e-6);
}

/**
 * Test the calculation of the generalized madelung matrix
 */
TEST(MadelungTests, CreateLattice1) {
  Matrix<double> vlat(3, 5);
  std::vector<double> vsq(5);

  std::vector<double> vn(3, 0.0);

  vsq[0] = 1.0;
  vsq[1] = 2.0;
  vsq[2] = 3.0;

  auto vnsq = 4.0;
  auto nv = 3;

  lsms::insert_ordered(vlat, vsq, nv, vn, vnsq);

  EXPECT_EQ(vsq[0], 1.0);
  EXPECT_EQ(vsq[1], 2.0);
  EXPECT_EQ(vsq[2], 3.0);
  EXPECT_EQ(vsq[3], 4.0);
}

TEST(MadelungTests, CreateLattice2) {
  Matrix<double> vlat(3, 6);
  std::vector<double> vsq(6);

  std::vector<double> vn(3, 0.0);

  vsq[0] = 1.0;
  vsq[1] = 2.0;
  vsq[2] = 3.0;
  vsq[3] = 5.0;
  vsq[4] = 6.0;

  auto vnsq = 4.0;
  auto nv = 5;

  lsms::insert_ordered(vlat, vsq, nv, vn, vnsq);

  EXPECT_EQ(vsq[0], 1.0);
  EXPECT_EQ(vsq[1], 2.0);
  EXPECT_EQ(vsq[2], 3.0);
  EXPECT_EQ(vsq[3], 4.0);
  EXPECT_EQ(vsq[4], 5.0);
  EXPECT_EQ(vsq[5], 6.0);
}

