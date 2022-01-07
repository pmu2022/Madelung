//
// Created by F.Moitzi on 16.12.2021.
//

#include "MultipoleMadelung.hpp"

#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>
#include <vector>
#include <complex>

#include "Matrix.hpp"
#include "Array3d.hpp"


#include "madelung.hpp"
#include "lattice_utils.hpp"

lsms::MultipoleMadelung::MultipoleMadelung(Matrix<double> lattice,
                                           Matrix<double> atom_position,
                                           int lmax,
                                           std::vector<int> global_position_index) :
    lattice{lattice},
    atom_position{atom_position},
    lmax{lmax},
    global_position_index{global_position_index} {

  num_atoms = atom_position.n_col(); // NOLINT(cppcoreguidelines-narrowing-conversions)
  local_num_atoms = global_position_index.size(); // NOLINT(cppcoreguidelines-narrowing-conversions)

  int jmax = (lmax + 1) * (lmax + 2) / 2;
  int kmax = (lmax + 1) * (lmax + 1);

  auto r_brav = lattice;
  auto k_brav = lattice;

  // 1. Scaling factors and rescale lattices and atomic positions
  auto scaling_factor = lsms::scaling_factor(lattice, lmax);

  r_brav.scale(1.0 / scaling_factor);
  atom_position.scale(1.0 / scaling_factor);
  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  auto omegbra = lsms::omega(r_brav);
  auto alat = scaling_factor * std::cbrt(3.0 * omegbra / (4.0 * M_PI * num_atoms));

  // 2. Calculate truncation spheres
  std::vector<int> r_nm(3);
  std::vector<int> k_nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  auto eta = lsms::calculate_eta(r_brav);


  lsms::real_space_trunc(r_brav, lmax, eta, rscut, r_nm);

  auto nrslat = num_latt_vectors(r_brav, rscut, r_nm);

  reciprocal_space_trunc(k_brav, lmax, eta, kncut, k_nm);
  auto nknlat = num_latt_vectors(k_brav, kncut, k_nm);

  // 3. Create the lattices
  Matrix<double> rslat;
  std::vector<double> rslatsq;

  Matrix<double> knlat;
  std::vector<double> knlatsq;

  std::tie(rslat, rslatsq) = lsms::create_lattice_and_sq(r_brav, rscut, r_nm, nrslat);
  std::tie(knlat, knlatsq) = lsms::create_lattice_and_sq(k_brav, kncut, k_nm, nknlat);

  auto omega = lsms::omega(r_brav);

  // 4. Calculate the madelung matrix
  madsum = Matrix<double>(num_atoms, local_num_atoms);

  // 5. Calculate DL matrix
  dl_matrix = Array3d<std::complex<double>>(num_atoms, kmax, local_num_atoms);

  for (auto i = 0; i < num_atoms; i++) {

    lsms::calculate_madelung(
        madsum,
        dl_matrix,
        atom_position,
        num_atoms,
        num_atoms,
        i, // id
        i, // myid
        jmax,
        kmax,
        lmax,
        omega,
        eta,
        scaling_factor,
        alat,
        nrslat,
        rslat,
        rslatsq,
        nknlat,
        knlat,
        knlatsq);

  }

  // 6. Dl factors
  dl_factor = lsms::calculate_dl_factor(kmax, jmax, lmax);

}

double lsms::MultipoleMadelung::getMadSum(int i, int j) const{
  return madsum(i, j);
}

std::complex<double> lsms::MultipoleMadelung::getDlMatrix(int i, int k, int j) {
  return dl_matrix(i, k, j);
}

double lsms::MultipoleMadelung::getDlFactor(int i, int j) {
  return dl_factor(i, j);
}




