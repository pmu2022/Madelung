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
 * Tests for the development of the Madelung potential
 */
TEST(MadelungTests, BasicLegacy) {
  int num_local_atoms = 2;
  int num_atoms = 2;
  int gindex[] = {1, 2};
  int lmax_rho = 3;
  int lmax_pot = 3;
  int iprint = 3;

  double bravais[] = {
      1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
  };

  double position[] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5};

  initMadelung(num_local_atoms, num_atoms, gindex, lmax_rho, lmax_pot, bravais,
               position, 10);

  printMadelungMatrix(&iprint);
  endMadelung();
}

/**
 * Tests for the development of the Madelung potential
 */
TEST(MadelungTests, BasicLegacy2) {
  int num_local_atoms = 2;
  int num_atoms = 2;
  int gindex[] = {1, 2};
  int lmax_rho = 3;
  int lmax_pot = 3;
  int iprint = 3;

  double bravais[] = {
      1.1, 0.2, 0.5, -0.1, 1.1, 0.1, 0.0, 0.0, 1.6,
  };

  double position[] = {0.0, 0.0, 0.1, 0.4, 0.4, 0.6};

  initMadelung(num_local_atoms, num_atoms, gindex, lmax_rho, lmax_pot, bravais,
               position, 10);

  printMadelungMatrix(&iprint);
  endMadelung();
}

/**
 * Check real space truncuation sphere
 */
TEST(MadelungTests, RealSpaceTrunc1) {
  int lmax = 3;
  Matrix<double> bravais(3, 3);
  bravais = 0.0;
  bravais(0, 0) = 3.5822399880299032;
  bravais(1, 1) = 3.5822399880299032;
  bravais(2, 2) = 3.5822399880299032;

  double eta = 0.60;
  double rscut = 0.0;

  std::vector<int> nm(3, 0);

  real_space_trunc(bravais, lmax, eta, rscut, nm);

  EXPECT_NEAR(rscut, 12.409243328345438, 1e-6);
}

/**
 * Check number of lattice vectors
 */
TEST(MadelungTests, NumberOfLatticeVectors) {
  Matrix<double> bravais(3, 3);
  bravais = 0.0;

  bravais(0, 0) = 6.28319;
  bravais(1, 1) = 6.28319;
  bravais(2, 2) = 6.28319;

  std::vector<int> nm(3, 8);

  double rscut = 48.75;

  auto nrslat_diff = num_latt_vectors(bravais, rscut, nm);
  EXPECT_EQ(nrslat_diff, 3145);
}

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

/**
 *
 */
TEST(MadelungTests, RealSpaceAndReciprocalSpace) {
  int lmax = 3;
  auto eta = 0.6;

  Matrix<double> bravais(3, 3);
  bravais = 0.0;

  bravais(0, 0) = 1.0;
  bravais(1, 1) = 1.0;
  bravais(2, 2) = 1.0;

  auto scaling_factor = lsms::scaling_factor(bravais, lmax);
  EXPECT_NEAR(scaling_factor, 0.81915494309189596, 1e-6);

  auto r_brav = bravais;
  auto k_brav = bravais;
  r_brav.scale(1.0 / scaling_factor);

  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  std::vector<int> nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  real_space_trunc(r_brav, lmax, eta, rscut, nm);
  auto nrslat = num_latt_vectors(r_brav, rscut, nm);

  reciprocal_space_trunc(k_brav, lmax, eta, kncut, nm);
  auto nknlat = num_latt_vectors(k_brav, kncut, nm);

  EXPECT_EQ(nrslat, 729);
  EXPECT_EQ(nknlat, 729);

  // Create lattice vectors for real space grid
  auto rslats = lsms::create_lattice(r_brav, rscut, nm, nrslat);

  EXPECT_NEAR(rslats(0, 0), 0.0, 1.0e-14);
  EXPECT_NEAR(rslats(1, 0), 0.0, 1.0e-14);
  EXPECT_NEAR(rslats(2, 0), 0.0, 1.0e-14);

  // Create lattice vectors for reciprocal space grid
  auto knlats = lsms::create_lattice(k_brav, kncut, nm, nknlat);
}

/**
 * Test the Madelung summations and DL matrix factors for l = 0
 *
 * This test is for a conventional body-centered cell
 *
 */
TEST(MadelungTests, MadelungSummations1) {
  // Global and local is the same
  int num_atoms = 2;
  int lmax = 3;
  int jmax = (lmax + 1) * (lmax + 2) / 2;
  int kmax = (lmax + 1) * (lmax + 1);
  auto eta = 0.6;

  Matrix<double> bravais(3, 3);
  bravais = 0.0;

  bravais(0, 0) = 1.0;
  bravais(1, 1) = 1.0;
  bravais(2, 2) = 1.0;

  Matrix<double> atom_position(3, 2);
  atom_position = 0.0;
  atom_position(0, 1) = 0.5;
  atom_position(1, 1) = 0.5;
  atom_position(2, 1) = 0.5;

  // 1. Scaling factors and rescale lattices and atomic positions
  auto scaling_factor = lsms::scaling_factor(bravais, lmax);
  auto r_brav = bravais;
  auto k_brav = bravais;
  r_brav.scale(1.0 / scaling_factor);
  atom_position.scale(1.0 / scaling_factor);
  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  auto omegbra = lsms::omega(r_brav);
  auto alat =
      scaling_factor * std::cbrt(3.0 * omegbra / (4.0 * M_PI * num_atoms));

  // 2. Calculate truncation spheres
  std::vector<int> nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  real_space_trunc(r_brav, lmax, eta, rscut, nm);
  auto nrslat = num_latt_vectors(r_brav, rscut, nm);

  reciprocal_space_trunc(k_brav, lmax, eta, kncut, nm);
  auto nknlat = num_latt_vectors(k_brav, kncut, nm);

  // 3. Create the lattices
  Matrix<double> rslat;
  std::vector<double> rslatsq;

  Matrix<double> knlat;
  std::vector<double> knlatsq;

  std::tie(rslat, rslatsq) =
      lsms::create_lattice_and_sq(r_brav, rscut, nm, nrslat);
  std::tie(knlat, knlatsq) =
      lsms::create_lattice_and_sq(k_brav, kncut, nm, nknlat);

  auto omega = lsms::omega(r_brav);

  // 4. Calculate the madelung matrix
  Matrix<double> madsum(num_atoms, num_atoms);

  // 5. Calculate DL matrix
  Array3d<std::complex<double>> DL_matrix(num_atoms, kmax, num_atoms);

  for (auto i = 0; i < num_atoms; i++) {
    lsms::calculate_madelung(madsum, DL_matrix, atom_position, num_atoms,
                             num_atoms,
                             i,  // id
                             i,  // myid
                             jmax, kmax, lmax, omega, eta, scaling_factor, alat,
                             nrslat, rslat, rslatsq, nknlat, knlat, knlatsq);
  }

  EXPECT_NEAR(madsum(0, 0), -0.2837297479481e1, 1e-12);
  EXPECT_NEAR(madsum(1, 0), -0.8019359700280, 1e-12);
  EXPECT_NEAR(madsum(0, 1), -0.8019359700280, 1e-12);
  EXPECT_NEAR(madsum(1, 1), -0.2837297479481e1, 1e-12);

  EXPECT_NEAR(std::real(DL_matrix(0, 0, 0)), -10.057957687339862, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(1, 0, 0)), -2.8427889965116413, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(1, 0, 1)), -10.057957687339862, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(0, 0, 1)), -2.8427889965116413, 1e-12);

  EXPECT_NEAR(std::imag(DL_matrix(0, 0, 0)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 0, 0)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 0, 1)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(0, 0, 1)), 0.0, 1e-12);
}

/**
 * Test the Madelung summations and DL matrix factors for l = 0
 *
 * This test is for a conventional body-centered cell
 *
 */
TEST(MadelungTests, MadelungSummations2) {
  // Global and local is the same
  int num_atoms = 2;
  int lmax = 3;
  int jmax = (lmax + 1) * (lmax + 2) / 2;
  int kmax = (lmax + 1) * (lmax + 1);

  Matrix<double> bravais(3, 3);
  bravais = 0.0;

  bravais(0, 0) = 1.1;
  bravais(1, 0) = 0.2;
  bravais(2, 0) = 0.5;

  bravais(0, 1) = -0.1;
  bravais(1, 1) = 1.1;
  bravais(2, 1) = 0.1;

  bravais(0, 2) = 0.0;
  bravais(1, 2) = 0.0;
  bravais(2, 2) = 1.6;

  Matrix<double> atom_position(3, 2);
  atom_position = 0.0;

  atom_position(0, 0) = 0.0;
  atom_position(1, 0) = 0.0;
  atom_position(2, 0) = 0.1;

  atom_position(0, 1) = 0.4;
  atom_position(1, 1) = 0.4;
  atom_position(2, 1) = 0.6;

  auto r_brav = bravais;
  auto k_brav = bravais;

  // 1. Scaling factors and rescale lattices and atomic positions
  auto scaling_factor = lsms::scaling_factor(bravais, lmax);
  EXPECT_NEAR(0.93651137065361856, scaling_factor, 1e-12);

  r_brav.scale(1.0 / scaling_factor);
  atom_position.scale(1.0 / scaling_factor);
  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  // 1.1 Test bravais lattices
  EXPECT_NEAR(r_brav(0, 0), 1.1745719640673216, 1e-12);
  EXPECT_NEAR(r_brav(1, 0), 0.21355853892133117, 1e-12);
  EXPECT_NEAR(r_brav(2, 0), 0.53389634730332791, 1e-12);

  EXPECT_NEAR(k_brav(0, 0), 5.2623592947212741, 1e-12);
  EXPECT_NEAR(k_brav(1, 1), 5.2623592947212741, 1e-12);
  EXPECT_NEAR(k_brav(2, 2), 3.6776715525608914, 1e-12);

  auto omegbra = lsms::omega(r_brav);
  auto alat =
      scaling_factor * std::cbrt(3.0 * omegbra / (4.0 * M_PI * num_atoms));

  // 2. Calculate truncation spheres
  std::vector<int> r_nm(3);
  std::vector<int> k_nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  auto eta = lsms::calculate_eta(r_brav);
  real_space_trunc(r_brav, lmax, eta, rscut, r_nm);
  EXPECT_EQ(r_nm[0], 4);
  EXPECT_EQ(r_nm[1], 4);
  EXPECT_EQ(r_nm[2], 3);
  EXPECT_NEAR(rscut, 10.400980527028162, 1e-12);
  auto nrslat = num_latt_vectors(r_brav, rscut, r_nm);

  reciprocal_space_trunc(k_brav, lmax, eta, kncut, k_nm);
  EXPECT_EQ(k_nm[0], 4);
  EXPECT_EQ(k_nm[1], 4);
  EXPECT_EQ(k_nm[2], 5);
  EXPECT_NEAR(kncut, 41.162263498561067, 1e-12);
  auto nknlat = num_latt_vectors(k_brav, kncut, k_nm);

  // 3. Create the lattices
  Matrix<double> rslat;
  std::vector<double> rslatsq;

  Matrix<double> knlat;
  std::vector<double> knlatsq;

  EXPECT_EQ(567, nrslat);
  EXPECT_EQ(891, nknlat);

  std::tie(rslat, rslatsq) =
      lsms::create_lattice_and_sq(r_brav, rscut, r_nm, nrslat);
  std::tie(knlat, knlatsq) =
      lsms::create_lattice_and_sq(k_brav, kncut, k_nm, nknlat);

  auto omega = lsms::omega(r_brav);

  // 4. Calculate the madelung matrix
  Matrix<double> madsum(num_atoms, num_atoms);

  // 5. Calculate DL matrix
  Array3d<std::complex<double>> DL_matrix(num_atoms, kmax, num_atoms);

  for (auto i = 0; i < num_atoms; i++) {
    lsms::calculate_madelung(madsum, DL_matrix, atom_position, num_atoms,
                             num_atoms,
                             i,  // id
                             i,  // myid
                             jmax, kmax, lmax, omega, eta, scaling_factor, alat,
                             nrslat, rslat, rslatsq, nknlat, knlat, knlatsq);
  }

  // 6. Dl factors
  auto dl_factor = lsms::calculate_dl_factor(kmax, jmax, lmax);

  EXPECT_NEAR(madsum(0, 0), -2.2231877288887785, 1e-12);
  EXPECT_NEAR(madsum(1, 0), -0.30394933422241815, 1e-12);
  EXPECT_NEAR(madsum(0, 1), -0.30394933422241815, 1e-12);
  EXPECT_NEAR(madsum(1, 1), -2.2231877288887785, 1e-12);

  EXPECT_NEAR(std::real(DL_matrix(0, 0, 0)), -7.8809953027095991, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(1, 0, 0)), -1.0774723358453844, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(1, 0, 1)), -7.8809953027095991, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(0, 0, 1)), -1.0774723358453844, 1e-12);

  EXPECT_NEAR(std::imag(DL_matrix(0, 0, 0)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 0, 0)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 0, 1)), 0.0, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(0, 0, 1)), 0.0, 1e-12);

  /*
   * (0,0)
   */
  EXPECT_NEAR(std::real(DL_matrix(0, 4, 0)), -2.3714279877125302, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(0, 4, 0)), 0.36116299470502256, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(0, 5, 0)), 2.3614205023378783, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(0, 5, 0)), -1.2316304346268181, 1e-12);

  /*
   * (1,0)
   */

  EXPECT_NEAR(std::real(DL_matrix(1, 1, 0)), 1.3254798267045567, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 1, 0)), 1.1053955921413914, 1e-12);

  /*
   * (0,1)
   */
  EXPECT_NEAR(std::real(DL_matrix(0, 1, 1)), -1.3254798267045567, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(0, 1, 1)), -1.1053955921413914, 1e-12);

  /*
   * (1,1)
   */
  EXPECT_NEAR(std::real(DL_matrix(1, 4, 1)), -2.3714279877125302, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 4, 1)), 0.36116299470502256, 1e-12);
  EXPECT_NEAR(std::real(DL_matrix(1, 5, 1)), 2.3614205023378783, 1e-12);
  EXPECT_NEAR(std::imag(DL_matrix(1, 5, 1)), -1.2316304346268181, 1e-12);

  /*
   * DL factors
   */

  EXPECT_NEAR(0.28209479177387842, dl_factor(0, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959454E-002, dl_factor(1, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959579E-002, dl_factor(2, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959454E-002, dl_factor(3, 0), 1e-12);
  EXPECT_NEAR(9.4031597257959593E-002, dl_factor(0, 1), 1e-12);
  EXPECT_NEAR(2.4278854013157353E-002, dl_factor(1, 1), 1e-12);
  EXPECT_NEAR(2.8034805800224084E-002, dl_factor(2, 1), 1e-12);
  EXPECT_NEAR(2.4278854013157353E-002, dl_factor(3, 1), 1e-12);
}