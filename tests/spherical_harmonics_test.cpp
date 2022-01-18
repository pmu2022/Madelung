
#include <gtest/gtest.h>

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>

#include "spherical_harmonics.hpp"
#include "integer_factors.hpp"


extern "C" {
#include <complex.h>
#include "gsl/gsl_sf.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
}

#include "utils.hpp"

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
 *
 */
TEST(SphericalHarmonicsTest, Clm3) {


  int lmax = 6;
  auto *clm_local = new double[(lmax + 1) * (lmax + 2) / 2];

  calc_clm(&lmax, clm_local);

  std::vector<double> ref{
      0.28209479177387831,
      0.48860251190292020,
      -0.34549414947133572,
      0.63078313050504031,
      -0.25751613468212653,
      0.12875806734106326,
      0.74635266518023113,
      -0.21545345607610053,
      6.8132365095552191E-002,
      -2.7814921575518955E-002,
      0.84628437532163492,
      -0.18923493915151210,
      4.4603102903819303E-002,
      -1.1920680675222410E-002,
      4.2145970709045986E-003,
      0.93560257962738924,
      -0.17081687924064815,
      3.2281355871636191E-002,
      -6.5894041742255317E-003,
      1.5531374585246052E-003,
      -4.9114518882630519E-004,
      1.0171072362820552,
      -0.15694305382900609,
      2.4814875652103469E-002,
      -4.1358126086839114E-003,
      7.5509261979682157E-004,
      -1.6098628745551694E-004,
      4.6472738199140594E-005
  };

  for (int i = 0; i < ref.size(); i++) {
    EXPECT_NEAR(ref[i], clm_local[i], 1.0e-12);
  }


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

/**
 * Compare the to analytic results
 */
TEST(SphericalHarmonics, SphHarmTest4) {

  int lmax = 3;

  std::vector<std::complex<double>> ylm((lmax + 1) * (lmax + 1));
  std::vector<std::complex<double>> ylm_gsl((lmax + 1) * (lmax + 1));
  std::vector<double> vec(3);

  vec[0] = 1.0;
  vec[1] = 1.0;
  vec[2] = 1.0;

  auto x = vec[0];
  auto y = vec[1];
  auto z = vec[2];

  // 1. Calculate data
  sph_harm_1(vec.data(), &lmax, ylm.data());


  auto r = norm(std::begin(vec), std::end(vec));
  auto theta = std::atan2(vec[1], vec[0]);
  auto cos_phi = vec[2] / r; //

  using namespace std::complex_literals;
  using namespace lsms;


  // reference
  auto yref_l2mm1 = 0.5 * std::sqrt(15.0 / (M_PI * 2.0)) * (x - y * 1i) * z / (r * r);
  auto yref_l2m1 = -0.5 * std::sqrt(15.0 / (M_PI * 2.0)) * (x + y * 1i) * z / (r * r);

  EXPECT_NEAR(std::real(yref_l2mm1), std::real(ylm[5]), 1e-12);
  EXPECT_NEAR(std::imag(yref_l2mm1), std::imag(ylm[5]), 1e-12);

  EXPECT_NEAR(std::real(yref_l2m1), std::real(ylm[7]), 1e-12);
  EXPECT_NEAR(std::imag(yref_l2m1), std::imag(ylm[7]), 1e-12);

}




/**
 * Test spherical harmonics that are used for the calculation of the Madelung constant
 * and reduced Green's function
 *
 * This is a comparison tests against the GSL implementation
 */

namespace lsms {

  gsl_complex convert_complex(std::complex<double> &z) {
    gsl_complex zr;
    GSL_SET_COMPLEX(&zr, std::real(z), std::imag(z));
    return zr;
  }

}

TEST(SphericalHarmonics, SphHarmTest3) {

  int lmax = 3;

  std::vector<std::complex<double>> ylm((lmax + 1) * (lmax + 1));
  std::vector<std::complex<double>> ylm_gsl((lmax + 1) * (lmax + 1));
  std::vector<double> vec(3);

  vec[0] = 1.0;
  vec[1] = 1.0;
  vec[2] = 1.0;

  // 1. Calculate data
  sph_harm_1(vec.data(), &lmax, ylm.data());


  auto r = norm(std::begin(vec), std::end(vec));
  auto theta = std::atan2(vec[1], vec[0]);
  auto cos_phi = vec[2] / r; //

  using namespace std::complex_literals;
  using namespace lsms;

  auto legendresize = gsl_sf_legendre_array_n(lmax);
  std::vector<double> ylm_compare(legendresize);

  auto legendre = gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SPHARM,
                                          lmax,
                                          cos_phi,
                                          -1.0,
                                          ylm_compare.data());

  auto lofk = lsms::get_lofk(lmax);
  auto mofk = lsms::get_mofk(lmax);

  auto kmax = lsms::get_kmax(lmax);

  for (auto k = 0; k < kmax; k++) {

    auto m = mofk[k];
    auto l = lofk[k];

    auto m_abs = std::abs(m);

    auto index = gsl_sf_legendre_array_index(l, m_abs);
    std::complex<double> z = m * theta * 1i;
    auto phase_factor = gsl_complex_exp(convert_complex(z));
    gsl_complex zr;
    GSL_SET_COMPLEX(&zr, std::real(ylm_compare[index]), 0.0);
    auto res = gsl_complex_mul(phase_factor, zr);

    ylm_gsl[k] = std::complex<double>(res.dat[0], res.dat[1]);

  }


  std::size_t i = 0;
  for (auto l = 0; l <= lmax; l++) {
    for (auto m = -l; m <= l; m++) {
      std::cout << l << " " << m << std::endl;
      std::cout << ylm_gsl[i] << std::endl;
      std::cout << ylm[i] << std::endl;
      i++;
    }
  }


}

