
#include <gtest/gtest.h>

#include <cmath>

#include "Array3d.hpp"
#include "Matrix.hpp"
#include "gaunt_factor.hpp"
#include "gaussq.h"

/**
 * Test the Gauss Quadratur
 */
TEST(GaussqTest, GaussLegendreQuadrature) {

  int kind = 1;
  int n = 2;
  int kpts = 0;
  double x1 = -1.0;
  double x2 = 1.0;

  auto *t = new double[n];
  auto *w = new double[n];

  gaussq(&kind, &n, &kpts, &x1, &x2, t, w);

  EXPECT_NEAR(-1.0 / sqrt(3.0), t[0], 1.0e-12);
  EXPECT_NEAR(1.0 / sqrt(3.0), t[1], 1.0e-12);

  EXPECT_NEAR(1.0, w[0], 1.0e-12);
  EXPECT_NEAR(1.0, w[1], 1.0e-12);

}


/**
 * Test the Gaunt factors
 */
TEST(GauntFactorTest, Lmax3) {


  /*
   * Reference results:
   *
   * MaxJ3 = lmax + 1
   *
   * cgnt(1:MaxJ3, 1:kmax, 1:kmax)
   *
   * cgnt(1, 1, 1)
   * cgnt(1, 2, 1)
   * cgnt(1, 1, 2)
   * cgnt(1, 2, 2)
   * cgnt(2, 2, 2)
   * cgnt(3, 3, 3)
   *
   *   0.28209479177387842
   *  -0.28209479177387836
   *   0.28209479177387836
   *  -0.12615662610100786
   *   0.28209479177387836
   *   0.0000000000000000
   *
   */

  int lmax = 3;
  int kmax = (lmax + 1) * (lmax + 1);

  Matrix<int> nj3(kmax, kmax);
  Array3d<int> kj3(lmax + 1, kmax, kmax);
  Array3d<double> cgnt(lmax + 1, kmax, kmax);

  gaunt_factor(&lmax, &cgnt[0], &kj3[0], &nj3[0]);

  EXPECT_NEAR(0.28209479177387842, cgnt(0, 0, 0), 1e-12);
  EXPECT_NEAR(-0.28209479177387842, cgnt(0, 1, 0), 1e-12);
  EXPECT_NEAR(0.28209479177387842, cgnt(0, 0, 1), 1e-12);
  EXPECT_NEAR(-0.12615662610100786, cgnt(0, 1, 1), 1e-12);
  EXPECT_NEAR(0.28209479177387842, cgnt(1, 1, 1), 1e-12);
  EXPECT_NEAR(0.0000000000000000, cgnt(2, 2, 2), 1e-12);

}


