//
// Created by F.Moitzi on 16.12.2021.
//

#include "MultipoleMadelung.hpp"

#define _USE_MATH_DEFINES

#include "math.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>

#include "Matrix.hpp"
#include "Array3d.hpp"

#include "Lattice.hpp"
#include "common.hpp"
#include "utils.hpp"
#include "spherical_harmonics.hpp"
#include "integer_factors.hpp"
#include "gaunt_factor.hpp"


double lsms::scaling_factor(Matrix<double> &bravais, int lmax) {

  constexpr auto max_iter = 100;
  constexpr auto fstep = 0.02;

  double rscut = 0.0;
  double kncut = 0.0;

  // Scaled bravais lattice
  Matrix<double> r_brav(3, 3);
  // Scaled reciprocal lattice
  Matrix<double> k_brav(3, 3);

  // Get the shortest axis
  auto a0 = sqrt(bravais(0, 0) * bravais(0, 0) +
                 bravais(1, 0) * bravais(1, 0) +
                 bravais(2, 0) * bravais(2, 0));

  auto a1 = sqrt(bravais(0, 1) * bravais(0, 1) +
                 bravais(1, 1) * bravais(1, 1) +
                 bravais(2, 1) * bravais(2, 1));

  auto a2 = sqrt(bravais(0, 2) * bravais(0, 2) +
                 bravais(1, 2) * bravais(1, 2) +
                 bravais(2, 2) * bravais(2, 2));


  // 1.0
  double scaling_fac = std::min({a0, a1, a2});

  auto eta = 0.5 + 0.1 * std::max({a0, a1, a2}) / scaling_fac;

  // 0.1591
  scaling_fac /= 2.0 * M_PI;


  for (int i = 0; i <= max_iter; i++) {

    r_brav = bravais;
    r_brav.scale(1 / scaling_fac);
    k_brav = 0.0;

    reciprocal_lattice(r_brav, k_brav, scaling_fac);

    std::vector<int> nm(3, 0);

    // Radius of real space truncation sphere
    lsms::real_space_trunc(r_brav, lmax, eta, rscut, nm);

    // Calculate number of lattice vectors
    auto nrslat = lsms::num_latt_vectors(r_brav, rscut, nm);

    // Radius of reciprocal space
    lsms::reciprocal_space_trunc(k_brav, lmax, eta, kncut, nm);
    // Calculate number of lattice vectors
    auto nknlat = lsms::num_latt_vectors(k_brav, kncut, nm);

    //printf("ALAT: %f %f %f\n", r_brav(0, 0), r_brav(1, 1), r_brav(2, 2));
    //printf("RS: %f KN: %f RSLAT: %d KNLAT: %d SC: %f\n", rscut, kncut, nrslat, nknlat, scaling_fac);


    if (nknlat < nrslat / 2) {
      scaling_fac = scaling_fac - fstep;
    } else if (nrslat < nknlat / 2) {
      scaling_fac = scaling_fac + fstep;
    } else {
      break;
    }


  }


  return scaling_fac;
}

void lsms::reciprocal_lattice(Matrix<double> &bravais,
                              Matrix<double> &reciprocal_bravais,
                              double &scale
) {

  // 248.050
  auto vol = (bravais(1, 0) * bravais(2, 1) -
              bravais(2, 0) * bravais(1, 1)) * bravais(0, 2) +
             (bravais(2, 0) * bravais(0, 1) -
              bravais(0, 0) * bravais(2, 1)) *
             bravais(1, 2) +
             (bravais(0, 0) *
              bravais(1, 1) -
              bravais(1, 0) *
              bravais(0, 1)) *
             bravais(2, 2);

  // 0.02533
  auto factor = 2.0 * M_PI / std::fabs(vol);

  reciprocal_bravais(0, 0) = factor * (bravais(1, 1) * bravais(2, 2) - bravais(2, 1) * bravais(1, 2));
  reciprocal_bravais(1, 0) = factor * (bravais(2, 1) * bravais(0, 2) - bravais(0, 1) * bravais(2, 2));
  reciprocal_bravais(2, 0) = factor * (bravais(0, 1) * bravais(1, 2) - bravais(1, 1) * bravais(0, 2));

  reciprocal_bravais(0, 1) = factor * (bravais(1, 2) * bravais(2, 0) - bravais(2, 2) * bravais(1, 0));
  reciprocal_bravais(1, 1) = factor * (bravais(2, 2) * bravais(0, 0) - bravais(0, 2) * bravais(2, 0));
  reciprocal_bravais(2, 1) = factor * (bravais(0, 2) * bravais(1, 0) - bravais(1, 2) * bravais(0, 0));

  reciprocal_bravais(0, 2) = factor * (bravais(1, 0) * bravais(2, 1) - bravais(2, 0) * bravais(1, 1));
  reciprocal_bravais(1, 2) = factor * (bravais(2, 0) * bravais(0, 1) - bravais(0, 0) * bravais(2, 1));
  reciprocal_bravais(2, 2) = factor * (bravais(0, 0) * bravais(1, 1) - bravais(1, 0) * bravais(0, 1));

}


void lsms::real_space_trunc(Matrix<double> &brav,
                            int lmax_mad,
                            double eta,
                            double &rscut,
                            std::vector<int> &nm) {

  rscut = 0.0;


  std::vector<double> r(3, 0.0);


  for (int i = 0; i < 3; i++) {
    r[i] = std::sqrt(
        brav(0, i) * brav(0, i) +
        brav(1, i) * brav(1, i) +
        brav(2, i) * brav(2, i));

  }


  for (int i = 0; i < 3; i++) {
    nm[i] = 0;

    auto term = 1.0;
    while (term > 0.5 * EPSI) {
      nm[i]++;
      auto gamma = gamma_func(nm[i] * r[i] / eta, lmax_mad);
      term = gamma[lmax_mad] / std::pow((nm[i] * r[i] / 2.0), lmax_mad + 1);
      // std::cout << term << " " << nm[i] << " " << eta << " " << lmax_mad << " " << r[i] << std::endl;
    }

  }


  for (int i = -1; i <= 1; i++) {

    for (int idx = 0; idx < 3; idx++) {
      r[idx] = i * brav(idx, 0) * nm[0];
    }

    for (int j = -1; j <= 1; j++) {

      for (int idx = 0; idx < 3; idx++) {
        r[idx] = r[idx] + j * brav(idx, 1) * nm[1];
      }

      for (int k = -1; k <= 1; k++) {

        for (int idx = 0; idx < 3; idx++) {
          r[idx] = r[idx] + k * brav(idx, 2) * nm[2];
        }

        rscut = std::max(rscut, norm(r.begin(), r.end()));

      }
    }
  }


}


void lsms::reciprocal_space_trunc(Matrix<double> &brav,
                                  int lmax_mad,
                                  double eta,
                                  double &kncut,
                                  std::vector<int> &nm) {

  auto fac = eta * eta / 4.0;

  kncut = 0.0;

  std::vector<double> r(3, 0.0);


  for (int i = 0; i < 3; i++) {
    r[i] = brav(0, i) * brav(0, i) +
           brav(1, i) * brav(1, i) +
           brav(2, i) * brav(2, i);

  }


  for (int i = 0; i < 3; i++) {
    nm[i] = 0;

    auto term = 1.0;
    while (term > 0.5 * EPSI) {
      nm[i]++;
      auto rm = nm[i] * nm[i] * r[i];
      term = std::exp(-fac * rm) * std::pow(std::sqrt(rm), lmax_mad - 2);
    }

  }


  for (int i = -1; i <= 1; i++) {

    for (int idx = 0; idx < 3; idx++) {
      r[idx] = i * brav(idx, 0) * nm[0];
    }

    for (int j = -1; j <= 1; j++) {

      for (int idx = 0; idx < 3; idx++) {
        r[idx] = r[idx] + j * brav(idx, 1) * nm[1];
      }

      for (int k = -1; k <= 1; k++) {

        for (int idx = 0; idx < 3; idx++) {
          r[idx] = r[idx] + k * brav(idx, 2) * nm[2];
        }

        kncut = std::max(kncut, norm(r.begin(), r.end()));

      }
    }
  }

}


int lsms::num_latt_vectors(Matrix<double> &brav, double cut, std::vector<int> &nm) {

  int number = 0;

  // To also include vectors on the boarder
  double vcut2 = cut * cut + 1e-6;

  std::vector<double> vn(3, 0.0);

  for (int x = -nm[0]; x <= nm[0]; x++) {
    for (int y = -nm[1]; y <= nm[1]; y++) {
      for (int z = -nm[2]; z <= nm[2]; z++) {

        vn[0] = x * brav(0, 0) + y * brav(0, 1) + z * brav(0, 2);
        vn[1] = x * brav(1, 0) + y * brav(1, 1) + z * brav(0, 2);
        vn[2] = x * brav(2, 0) + y * brav(2, 1) + z * brav(0, 2);

        auto norm = norm_sq(vn.begin(), vn.end());

        if (norm <= vcut2) {
          number++;
        }

      }
    }
  }

  return number;

}


__attribute__((unused)) int lsms::num_latt_vectors_legacy(Matrix<double> &brav, double cut, std::vector<int> &nm) {


  int number = 0;

  // To also include vectors on the boarder
  double vcut2 = cut * cut + 1e-6;

  int tn1p1 = 2 * nm[0] + 1;
  int tn2p1 = 2 * nm[1] + 1;
  int tn3p1 = 2 * nm[2] + 1;

  int nt12 = tn1p1 * tn2p1;
  int nt123 = nt12 * tn3p1;

  int n1, n2, n3;

  std::vector<double> vn(3, 0.0);

  for (int i = 0; i < nt123; i++) {

    n1 = i - 1;
    n1 = n1 % tn1p1 - nm[0];
    n2 = (i - 1) / tn1p1;
    n2 = n2 % tn2p1 - nm[1];
    n3 = (i - 1) / nt12;
    n3 = n3 % tn3p1 - nm[2];

    vn[0] = n1 * brav(0, 0) + n2 * brav(0, 1) + n3 * brav(0, 2);
    vn[1] = n1 * brav(1, 0) + n2 * brav(1, 1) + n3 * brav(1, 2);
    vn[2] = n1 * brav(2, 0) + n2 * brav(2, 1) + n3 * brav(2, 2);

    if (norm_sq(vn.begin(), vn.end()) <= vcut2) {
      number++;
    }

  }

  return number;


}


void lsms::generate_lattice(Matrix<double> &r_brav,
                            Matrix<double> &k_brav,
                            int lmax,
                            double eta) {

  std::vector<int> nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  Matrix<double> rslat;
  std::vector<double> rslat_sq;

  Matrix<double> knlat;
  std::vector<double> knlat_sq;

  // Radius of real space truncation sphere
  lsms::real_space_trunc(r_brav, lmax, eta, rscut, nm);
  // Calculate number of lattice vectors
  auto nrslat = lsms::num_latt_vectors(r_brav, rscut, nm);
  // Create lattice vectors for real space grid
  tie(rslat, rslat_sq) = lsms::create_lattice_and_sq(r_brav, rscut, nm, nrslat);

  // Radius of reciprocal space
  lsms::reciprocal_space_trunc(k_brav, lmax, eta, kncut, nm);
  // Calculate number of lattice vectors
  auto nknlat = lsms::num_latt_vectors(k_brav, kncut, nm);
  // Create lattice vectors for reciprocal space grid
  tie(knlat, knlat_sq) = lsms::create_lattice_and_sq(k_brav, kncut, nm, nknlat);


}


void lsms::init_madelung(Matrix<double> &bravais,
                         Matrix<double> &atom_position,
                         int lmax) {

  /*
   * Calculate the scale
   */

  double scale = scaling_factor(bravais, lmax);

  auto r_brav = bravais;
  auto k_brav = bravais;

  /*
   *  change units so that both vbrar and atom_posi_* are in
   *  in the units of scale
   */
  r_brav.scale(1.0 / scale);
  atom_position.scale(1.0 / scale);

  reciprocal_lattice(r_brav, k_brav, scale);


}


/**
 * Reciprocal space term of Madelung sum:
 *
 *                                       2   2       -> ->
 *                4*pi          1    -eta *Kq /4 - i*Kq*aij
 *       term1 =  ---- * sum  ----- e
 *                tau    q<>0    2
 *                             Kq
 *
 *       note sum starts at 2 since kn=0.0 of 1/(rs-aij) is
 *       canceled by kn=0.0 of the 1/rs sum.
 *
 */
static double reciprocal_space_term(Matrix<double> &knlat,
                                    std::vector<double> &knlatsq,
                                    std::vector<double> &aij,
                                    int nknlat,
                                    double eta,
                                    double omega) {

  auto term12 = 0.0;
  auto fac = eta * eta / 4.0;

  for (int i = nknlat - 1; i > 0; i--) {
    term12 += std::exp(-fac * knlatsq[i]) / knlatsq[i]
              * std::cos(knlat(0, i) * aij[0]
                         + knlat(1, i) * aij[1]
                         + knlat(2, i) * aij[2]);
  }
  term12 = 4.0 * M_PI / omega * term12;

  return term12;
}

/**
 * Real space term of Madelung sum:
 *
 *                              ->   ->
 *                     1 - erf(|Rn + aij|/eta)
 *        term2 = sum -------------------------
 *                 n         ->   ->
 *                         |Rn + aij|
 *
 *        note for calculation of aij=0.0 term ibegin=2.
 *
 */
template<class T>
static T real_space_term(Matrix<T> &rslat,
                         std::vector<T> &aij,
                         int nrslat,
                         int ibegin,
                         T eta) {

  /*
  *  subtract aij from rslat and calculate rslatmd which is used in
  *  calculating the real-space integral
  *  rslatmd, and aij are in the units of a0 = 1
  */

  std::vector<T> rslatmd(nrslat);

  for (auto i = 0; i < nrslat; i++) {

    rslatmd[i] = std::sqrt(
        (rslat(0, i) - aij[0]) * (rslat(0, i) - aij[0])
        + (rslat(1, i) - aij[1]) * (rslat(1, i) - aij[1])
        + (rslat(2, i) - aij[2]) * (rslat(2, i) - aij[2]));


  }

  auto rterm = 0.0;

  for (auto i = nrslat - 1; i >= ibegin; i--) {
    rterm += std::erfc(rslatmd[i] / eta) / rslatmd[i];
  }


  return rterm;

}


/**
 *
 * Dl sum
 *
 */
static std::vector<std::complex<double>> dlsum(
    std::vector<double> &aij,
    Matrix<double> &rslat,
    int nrslat,
    int ibegin,
    Matrix<double> &knlat,
    int nknlat,
    double omega,
    int lmax_mad,
    int kmax_mad,
    double eta) {

  std::vector<std::complex<double>> Ylm(kmax_mad, std::complex<double>(0.0, 0.0));
  std::vector<double> vec(3);
  std::vector<std::complex<double>> dlm(kmax_mad, 0.0);

  auto aij2 = std::inner_product(aij.begin(), aij.end(), aij.begin(), 0.0);
  auto lofk = lsms::get_lofk(lmax_mad);

  /*
   * The unit of DL_matrix needs to be resumed to a0 = alat later
   */
  for (int i = nrslat - 1; i >= ibegin; i--) {

    for (int j = 0; j < 3; j++) {
      vec[j] = aij[j] + rslat(j, i);
    }

    auto vlen = norm(vec.begin(), vec.end());

    // Ylm
    sph_harm_1(vec.data(), &lmax_mad, Ylm.data());

    // Gamma
    auto gamma_l = gamma_func(vlen / eta, lmax_mad);

    auto vhalf = 0.5 * vlen;

    for (auto kl = kmax_mad - 1; kl > 0; kl--) {
      auto l = lofk[kl];
      dlm[kl] = dlm[kl] + gamma_l[l] * Ylm[kl] / std::pow(vhalf, l + 1);
    }

    // std::printf("%16.12f %16.12f %16.12f %16.12f %16.12e\n", vec[0], vec[1], vec[2], vlen, gamma_l[lofk[kmax_mad - 1]]);
  }

  auto rfac = 4.0 * std::sqrt(M_PI);
  for (int i = 1; i < kmax_mad; ++i) {
    dlm[i] *= rfac;
  }

  /*
   *
   */
  rfac = -eta * eta / 4.0;

  std::vector<std::complex<double>> ctmp(kmax_mad, 0.0);

  for (int i = nknlat - 1; i > 0; i--) {

    for (int j = 0; j < 3; j++) {
      vec[j] = knlat(j, i);
    }

    auto vlen = norm(vec.begin(), vec.end());
    auto knlatsq = norm_sq(vec.begin(), vec.end());

    //std::printf("%16.12f %16.12f %16.12f %16.12f %16.12f\n", vec[0], vec[1], vec[2], vlen, knlatsq);

    // Ylm
    sph_harm_1(vec.data(), &lmax_mad, Ylm.data());

    auto expfac = std::exp(rfac * knlatsq) / knlatsq;


    double tfac;
    double sintfac, costfac;

    if (aij2 > 1e-8) {
      tfac = -(knlat(0, i) * aij[0] + knlat(1, i) * aij[1] + knlat(2, i) * aij[2]);
      sintfac = sin(tfac);
      costfac = cos(tfac);
    } else {
      tfac = 0.0;
      sintfac = 0.0;
      costfac = 1.0;
    }

    auto cfac = -vlen;
    for (auto l = 1; l <= lmax_mad; l += 2) {
      for (auto m = -l; m <= l; m++) {
        auto kl = (l + 1) * (l + 1) - l + m - 1;
        ctmp[kl] = ctmp[kl] + expfac * cfac * sintfac * Ylm[kl];
      }
      cfac = -cfac * knlatsq;
    }

    cfac = -knlatsq;
    for (auto l = 2; l <= lmax_mad; l += 2) {
      for (auto m = -l; m <= l; m++) {
        auto kl = (l + 1) * (l + 1) - l + m - 1;
        ctmp[kl] = ctmp[kl] + expfac * cfac * costfac * Ylm[kl];
      }
      cfac = -cfac * knlatsq;
    }

  }

  rfac = 16.0 * M_PI * M_PI / omega;

  for (auto i = 1; i < kmax_mad; i++) {
    dlm[i] = dlm[i] + rfac * ctmp[i];
  }

  return dlm;

}

/**
 * Calculate the madelung constants
 */
void lsms::calculate_madelung(
    Matrix<double> &madmat,
    Array3d<std::complex<double>> &DL_matrix,
    Matrix<double> &atom_position,
    int global_natoms,
    int local_natoms,
    int id,
    int myatom,
    int jmax_mad,
    int kmax_mad,
    int lmax_mad,
    double omega,
    double eta,
    double a0,
    double alat,
    int nrslat,
    Matrix<double> &rslat,
    std::vector<double> &rslatsq,
    int nknlat,
    Matrix<double> &knlat,
    std::vector<double> &knlatsq) {


  // Auxiliar variables
  std::vector<double> aij(3);
  double r0tm;
  int ibegin;


  // First terms
  auto term0 = -M_PI * eta * eta / omega;


  for (int i = 0; i < global_natoms; i++) {

    /*
     * a_ij in unit of a0
     */
    for (int idx = 0; idx < 3; idx++) {
      aij[idx] = atom_position(idx, myatom) - atom_position(idx, i);
    }

    // Real space terms: first terms
    if (i == myatom) {
      ibegin = 1;
      r0tm = -2.0 / std::sqrt(M_PI) / eta;
    } else {
      ibegin = 0;
      r0tm = 0.0;
    }

    // Reciprocal space term
    auto term1 = reciprocal_space_term(knlat,
                                       knlatsq,
                                       aij,
                                       nknlat,
                                       eta,
                                       omega);


    // Real space term
    auto term2 = real_space_term(rslat, aij, nrslat, ibegin, eta);

    madmat(i, id) = term1 + term2 + r0tm + term0;
    madmat(i, id) = madmat(i, id) / a0;


    if (jmax_mad > 1) {

      /*
       * DL_factors(global_natoms, k, local_natoms)
       */

      // 1. First factor for k = 0
      DL_matrix(i, 0, id) = madmat(i, id) * Y0inv;

      std::vector<std::complex<double>> dlm(kmax_mad, 0.0);

      auto lofk = lsms::get_lofk(kmax_mad);

      dlm = dlsum(aij, rslat, nrslat, ibegin, knlat, nknlat, omega, lmax_mad, kmax_mad, eta);

      // 2. Calculate all other factors
      for (int kl = 1; kl < kmax_mad; kl++) {

        auto l = lofk[kl];
        DL_matrix(i, kl, id) = dlm[kl] * std::pow(alat / a0, l) / a0;

      }


    }


  }


}

double lsms::calculate_eta(Matrix<double> &bravais) {

  auto a0 = sqrt(bravais(0, 0) * bravais(0, 0) +
                 bravais(1, 0) * bravais(1, 0) +
                 bravais(2, 0) * bravais(2, 0));

  auto a1 = sqrt(bravais(0, 1) * bravais(0, 1) +
                 bravais(1, 1) * bravais(1, 1) +
                 bravais(2, 1) * bravais(2, 1));

  auto a2 = sqrt(bravais(0, 2) * bravais(0, 2) +
                 bravais(1, 2) * bravais(1, 2) +
                 bravais(2, 2) * bravais(2, 2));


  auto scaling_fac = std::min({a0, a1, a2});

  return 0.5 + 0.1 * std::max({a0, a1, a2}) / scaling_fac;

}

Matrix<double> lsms::calculate_dl_factor(int kmax_mad, int jmax_mad, int lmax_mad) {

  // Variable
  Matrix<double> dl_factor(kmax_mad, jmax_mad);

  std::vector<double> factmat(lmax_mad + 1);
  factmat[0] = 1.0;

  for (int l = 1; l <= lmax_mad; ++l) {
    factmat[l] = factmat[l - 1] / (2.0 * l + 1.0);
  }

  Array3d<double> cgnt(lmax_mad + 1,
                       (lmax_mad + 1) * (lmax_mad + 1),
                       (lmax_mad + 1) * (lmax_mad + 1));

  Matrix<int> nj3(kmax_mad, kmax_mad);
  Array3d<int> kj3(lmax_mad + 1, kmax_mad, kmax_mad);

  gaunt_factor(&lmax_mad, &cgnt[0], &kj3[0], &nj3[0]);

  auto lofj = lsms::get_lofj(lmax_mad);
  auto kofj = lsms::get_kofj(lmax_mad);

  auto lofk = lsms::get_lofk(lmax_mad);

  auto mofk = lsms::get_mofk(lmax_mad);
  auto mofj = lsms::get_mofj(lmax_mad);

  for (int jl_pot = 0; jl_pot < jmax_mad; jl_pot++) {
    auto l_pot = lofj[jl_pot];
    auto kl_pot = kofj[jl_pot];

    for (int kl_rho = 0; kl_rho < kmax_mad; kl_rho++) {
      auto l_rho = lofk[kl_rho];
      auto l_sum = l_pot + l_rho;

      auto m_dif = mofk[kl_rho] - mofj[jl_pot];
      auto kl = (l_sum + 1) * (l_sum + 1) - l_sum + m_dif - 1;

      auto gaunt = get_gaunt_factor(cgnt, nj3, kj3, kl_pot, kl_rho, kl);
      dl_factor(kl_rho, jl_pot) = gaunt * factmat[l_pot] * factmat[l_rho];

    }


  }


  return dl_factor;
}

std::tuple<Matrix<double>, Array3d<std::complex<double>>, Matrix<double> >
lsms::get_madelung(int num_atoms, int lmax, Matrix<double> &bravais, Matrix<double> &atom_position) {


  int jmax = (lmax + 1) * (lmax + 2) / 2;
  int kmax = (lmax + 1) * (lmax + 1);


  auto r_brav = bravais;
  auto k_brav = bravais;

  // 1. Scaling factors and rescale lattices and atomic positions
  auto scaling_factor = lsms::scaling_factor(bravais, lmax);

  r_brav.scale(1.0 / scaling_factor);
  atom_position.scale(1.0 / scaling_factor);
  lsms::reciprocal_lattice(r_brav, k_brav, scaling_factor);


  auto omegbra = lsms::omega(r_brav);
  auto alat = scaling_factor * std::cbrt(3.0 * omegbra / (4.0 * M_PI * num_atoms));

  // 2. Calculate truncation spheres
  std::vector<int> r_nm(3);
  std::vector<int> k_nm(3);
  double rscut = 0.0;
  double kncut = 0.0;

  auto eta = lsms::calculate_eta(r_brav);
  lsms::real_space_trunc(r_brav, lmax, eta, rscut, r_nm);
  auto nrslat = lsms::num_latt_vectors(r_brav, rscut, r_nm);

  lsms::reciprocal_space_trunc(k_brav, lmax, eta, kncut, k_nm);
  auto nknlat = lsms::num_latt_vectors(k_brav, kncut, k_nm);

  // 3. Create the lattices
  Matrix<double> rslat;
  std::vector<double> rslatsq;

  Matrix<double> knlat;
  std::vector<double> knlatsq;


  std::tie(rslat, rslatsq) = lsms::create_lattice_and_sq(r_brav, rscut, r_nm, nrslat);
  std::tie(knlat, knlatsq) = lsms::create_lattice_and_sq(k_brav, kncut, k_nm, nknlat);

  auto omega = lsms::omega(r_brav);

  // 4. Calculate the madelung matrix
  Matrix<double> madsum(num_atoms, num_atoms);

  // 5. Calculate DL matrix
  Array3d<std::complex<double>> DL_matrix(num_atoms, kmax, num_atoms);

  for (auto i = 0; i < num_atoms; i++) {

    lsms::calculate_madelung(
        madsum,
        DL_matrix,
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
  auto dl_factor = lsms::calculate_dl_factor(kmax, jmax, lmax);

  return std::make_tuple(madsum, DL_matrix, dl_factor);

}
