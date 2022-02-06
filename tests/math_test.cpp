/**
 *
 * Tests for the legendre functionals
 *
 */

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "integer_factors.hpp"
#include "spherical_harmonics.hpp"


TEST(MathTest, SphericalHarmonicsVector1) {

  using namespace lsms::math;

  auto lmax = 3;

  std::array<double, 3> vec{{1.0, 1.0, 1.0}};


  auto ylm = spherical_harmonics(vec, 3);

  // l = 0  m = 0
  const auto real_ref00 = sqrt(1.0 / (4.0 * M_PI));
  const auto imag_ref00 = 0.0;

  EXPECT_NEAR(real_ref00, std::real(ylm[0]), 1.0e-12);
  EXPECT_NEAR(imag_ref00, std::imag(ylm[0]), 1.0e-12);


}

TEST(MathTest, SphericalHarmonicsVector2) {

  using namespace lsms::math;

  auto lmax = 3;

  std::array<double, 3> vec{{-4.8830810748721376, -4.8830810748721376, -4.8830810748721376}};


  auto ylm = spherical_harmonics(vec, 3);

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


}


TEST(MathTest, SphericalHarmonicsVector3) {

  using namespace lsms::math;

  auto lmax = 30;

  std::array<double, 3> vec{{1.0, 1.0, 1.0}};

  auto ylm = spherical_harmonics(vec, lmax);

  using namespace std::complex_literals;
  using namespace lsms;

  auto lofk = lsms::get_lofk(lmax);
  auto mofk = lsms::get_mofk(lmax);
  auto kmax = lsms::get_kmax(lmax);

  const std::vector<std::complex<double>> compare_res_l30{
      -4.34262048809611e-18 + 0.0016114564738045761i,
      -0.006241144086183253 + 0.006241144086183253i,
      -0.03217445746136457 - 2.7581642441053956e-17i,
      -0.06315628610585003 - 0.06315628610585024i,
      -1.9386605149095688e-16 - 0.197754615196241i,
      0.25027973459899605 - 0.25027973459899605i,
      0.5067141193910657 + 3.72327494637023e-16i,
      0.3916518217721521 + 0.39165182177215335i,
      -9.713785172606613e-16 + 0.3964951752214758i,
      -0.0304294314775328 + 0.030429431477532815i,
      0.32322578670044794 + 1.9791871254229436e-16i,
      0.304949898509267 + 0.30494989850926796i,
      -9.108325125330795e-17 + 0.16527804192055603i,
      0.17692709022075687 - 0.176927090220757i,
      0.39546073364747425 + 1.9371988865993734e-16i,
      0.06810319540693736 + 0.06810319540693754i,
      1.3329964287978602e-16 - 0.31099262479503575i,
      0.22601061085123725 - 0.2260106108512375i,
      -0.09450624134527816 - 3.472102978888264e-17i,
      -0.2627318257447544 - 0.26273182574475507i,
      3.325528324220362e-17 - 0.10861999807734685i,
      -0.2221077634394482 + 0.2221077634394481i,
      -0.24396200057995396 - 5.975345662476507e-17i,
      0.1540482476931519 + 0.15404824769315198i,
      -5.775517454110126e-17 + 0.31440452643876676i,
      0.09013690777911247 - 0.09013690777911244i,
      0.34255235088485253 + 4.1950564005153565e-17i,
      -0.04451964580493865 - 0.044519645804938654i,
      2.1445150726182942e-17 - 0.3502258894746451i,
      -0.021430138556263315 + 0.02143013855626331i,
      -0.35125454366556247 + 0i,
      0.02143013855626331 + 0.021430138556263308i,
      2.1445150726182942e-17 + 0.3502258894746451i,
      0.04451964580493865 - 0.044519645804938654i,
      0.3425523508848526 - 4.195056400515357e-17i,
      -0.09013690777911247 - 0.09013690777911244i,
      -5.775517454110126e-17 - 0.31440452643876676i,
      -0.1540482476931519 + 0.15404824769315198i,
      -0.24396200057995399 + 5.975345662476508e-17i,
      0.22210776343944827 + 0.22210776343944816i,
      3.325528324220362e-17 + 0.10861999807734683i,
      0.2627318257447544 - 0.26273182574475507i,
      -0.09450624134527814 + 3.4721029788882634e-17i,
      -0.2260106108512373 - 0.22601061085123755i,
      1.33299642879786e-16 + 0.3109926247950357i,
      -0.06810319540693738 + 0.06810319540693756i,
      0.39546073364747425 - 1.9371988865993734e-16i,
      -0.17692709022075684 - 0.17692709022075698i,
      -9.1083251253308e-17 - 0.16527804192055612i,
      -0.30494989850926696 + 0.30494989850926785i,
      0.3232257867004479 - 1.9791871254229431e-16i,
      0.03042943147753279 + 0.030429431477532805i,
      -9.713785172606613e-16 - 0.3964951752214758i,
      -0.39165182177215224 + 0.3916518217721535i,
      0.5067141193910659 - 3.723274946370232e-16i,
      -0.2502797345989962 - 0.2502797345989963i,
      -1.9386605149095688e-16 + 0.197754615196241i,
      0.06315628610585003 - 0.06315628610585024i,
      -0.032174457461364575 + 2.7581642441053962e-17i,
      0.006241144086183255 + 0.006241144086183255i,
      -4.34262048809611e-18 - 0.0016114564738045761i,
  };

  std::size_t i = 0;
  for (auto k = 900; k < kmax; k++) {
    EXPECT_NEAR(std::real(compare_res_l30[i]), std::real(ylm[k]), 1e-14);
    EXPECT_NEAR(std::imag(compare_res_l30[i]), std::imag(ylm[k]), 1e-14);
    i++;
  }

}