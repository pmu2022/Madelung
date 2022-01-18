/**
 *
 * This tests cover utily routines

 *
 */

#include <gtest/gtest.h>

#define _USE_MATH_DEFINES


#include <tuple>
#include <vector>

#include "common.hpp"
#include "lattice_utils.hpp"
#include "madelung.hpp"



/**
 * Test the scaling factors
 */
TEST(MadelungTests, ScalingFactor1) {
  int lmax = 3;

  lsms::matrix<double> bravais(3, 3);
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
  lsms::matrix<double> vlat(3, 5);
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
  lsms::matrix<double> vlat(3, 6);
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

