/**
 * Tests the creation of the multipole moments
 *
 */

#include <gtest/gtest.h>

#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "MultipoleMoments.hpp"
#include "NDArray.hpp"
#include "integer_factors.hpp"
#include "integrator.hpp"
#include "utils.hpp"

TEST(MultipoleMomentsTests, IntegrationsReal) {
  std::size_t N = 10001;

  std::vector<double> rmesh = lsms::math::linspace(0, 10, N);
  std::vector<double> func(N);

  for (auto i = 0; i < N; i++) {
    func[i] = rmesh[i] * rmesh[i] + 10;
  }

  auto res = lsms::simpson_nonuniform(rmesh, func, N);
  EXPECT_NEAR(res, 433.3333333, 1e-6);
}

TEST(MultipoleMomentsTests, IntegrationsComplex) {
  std::size_t N = 10001;

  std::vector<double> rmesh = lsms::math::linspace(0, 10, N);
  std::vector<std::complex<double>> func(N);

  using namespace std::complex_literals;

  for (auto ir = 0; ir < N; ir++) {
    func[ir] = rmesh[ir] * rmesh[ir] + 10 + (rmesh[ir] * rmesh[ir] + 10) * 1i;
  }

  auto res = lsms::simpson_nonuniform(rmesh, func, N);

  EXPECT_NEAR(std::imag(res), 433.3333333, 1e-6);
  EXPECT_NEAR(std::real(res), 433.3333333, 1e-6);
}

TEST(MultipoleMomentsTests, ConstructMomentFromDenTildaL) {
  std::ifstream infile("RhoTildaL");

  std::vector<double> r_mesh;
  lsms::NDArray<std::complex<double>, 2> rhoTildaL(1036, 5);

  // Read data
  std::string line;

  auto idx = 0;
  while (std::getline(infile, line)) {
    std::istringstream iss(line);

    if (line[0] == '#') {
      continue;
    }

    int ind, spec;
    iss >> ind >> spec;

    double rmesh;
    iss >> rmesh;
    r_mesh.emplace_back(rmesh);

    double real;
    double imag;
    for (auto i = 0; i < 5; i++) {
      iss >> real >> imag;
      rhoTildaL(idx, i) = std::complex<double>(real, imag);
    }
    idx++;
  }

  auto nr_ps = 990;

  auto lmax = 7;
  auto jmax = lsms::get_jmax(lmax);

  lsms::NDArray<std::complex<double>, 2> density_tilda_l(1036, jmax);
  density_tilda_l = 0;

  std::vector<std::size_t> indices{0, 10, 14, 21, 25};

  // Construct the density
  auto jdx = 0;
  for (auto idx : indices) {
    for (auto ir = 0; ir < 1036; ir++) {
      density_tilda_l(ir, idx) = rhoTildaL(ir, jdx);
    }
    jdx++;
  }

  auto q_lm = lsms::multipole_moms::construct_multi_mom(r_mesh, density_tilda_l,
                                                        nr_ps, lmax);

  using namespace std::complex_literals;

  std::vector<std::complex<double>> ref = {
      (8.959096192752E+01 + 0.000000000000E+00i),
      (1.428840665798E-01 + 8.929488936022E-23i),
      (8.538956207453E-02 - 3.825172365396E-16i),
      (3.707204598109E-04 + 8.497617823544E-27i),
      (-6.935544732297E-04 + 7.546391765593E-17i)};

  auto i = 0;
  for (auto jl : indices) {
    EXPECT_NEAR(std::real(q_lm[jl]), std::real(ref[i]), 1e-5);
    EXPECT_NEAR(std::imag(q_lm[jl]), std::imag(ref[i]), 1e-5);
    i++;
  }

}
