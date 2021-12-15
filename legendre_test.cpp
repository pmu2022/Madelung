
#include <gtest/gtest.h>

#include <math.h>

#include "legendre.h"

TEST(LegendreTest, Lmax1) {

  int lmax = 1;
  double x = 0.2;

  auto *plm = new double[(lmax + 1) * (lmax + 2) / 2];
  legendre(&lmax, &x, plm);


  // l = 0  m = 0
  EXPECT_NEAR(1.0, plm[0], 1.0e-12);
  // l = 1  m = 0
  EXPECT_NEAR(x, plm[1], 1.0e-12);
  // l = 1  m = 1
  EXPECT_NEAR(sqrt(1.0 - x * x), plm[2], 1.0e-12);

  delete[] plm;

}

TEST(LegendreTest, Lmax2) {

  int lmax = 2;
  double x = 0.2;

  auto *plm = new double[(lmax + 1) * (lmax + 2) / 2];
  legendre(&lmax, &x, plm);


  // l = 0  m = 0
  EXPECT_NEAR(1.0, plm[0], 1.0e-12);
  // l = 1  m = 0
  EXPECT_NEAR(x, plm[1], 1.0e-12);
  // l = 1  m = 1
  EXPECT_NEAR(sqrt(1.0 - x * x), plm[2], 1.0e-12);
  // l = 2  m = 0
  EXPECT_NEAR(0.5 * ( 3.0 * x * x - 1.0), plm[3], 1.0e-12);
  // l = 2  m = 1
  EXPECT_NEAR(3 * x * sqrt(1.0 - x * x), plm[4], 1.0e-12);
  // l = 2  m = 2
  EXPECT_NEAR(3 * ( 1 - x * x), plm[5], 1.0e-12);


  delete[] plm;

}


