/**
 *
 * Tests for the legendre functionals
 *
 */

#include <gtest/gtest.h>
#include <math.h>

extern "C" {
#include "gsl/gsl_sf.h"
}

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
  EXPECT_NEAR(0.5 * (3.0 * x * x - 1.0), plm[3], 1.0e-12);
  // l = 2  m = 1
  EXPECT_NEAR(3 * x * sqrt(1.0 - x * x), plm[4], 1.0e-12);
  // l = 2  m = 2
  EXPECT_NEAR(3 * (1 - x * x), plm[5], 1.0e-12);

  delete[] plm;
}

TEST(LegendreTest, Lmax20) {
  int lmax = 20;

  for (double x : {0.00001, 0.001, 0.01, 0.1, 0.5, 0.95}) {
    std::vector<double> plm((lmax + 1) * (lmax + 2) / 2);
    legendre(&lmax, &x, plm.data());

    auto legendresize = gsl_sf_legendre_array_n(lmax);
    std::vector<double> plm_compare(legendresize);

    auto res = gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_NONE, lmax, x, 1.0,
                                       plm_compare.data());

    std::size_t i = 0;
    for (auto l = 0; l <= 20; l++) {
      for (auto m = 0; m <= l; m++) {
        auto index = gsl_sf_legendre_array_index(l, m);

        auto rel_error = std::abs(plm_compare[index] - plm[i]) /
                         std::abs(plm_compare[index]);
        EXPECT_TRUE(rel_error <= 1.0e-12);
        // std::cout << plm_compare[index] << " " << plm[i] << std::endl;
        i++;
      }
    }
  }
}
