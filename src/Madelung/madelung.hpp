//
// Created by F.Moitzi on 06.01.2022.
//

#ifndef SRC_MADELUNG_MADELUNG_HPP
#define SRC_MADELUNG_MADELUNG_HPP

#include <vector>
#include <complex>

#include "Matrix.hpp"
#include "Array3d.hpp"

namespace lsms {

  constexpr auto EPSI = 1e-14;

  /**
   * Calculate the scaling factor to get a balanced number
   * of reciprocal and real-space lattice vectors
   */
  double scaling_factor(const Matrix<double> &bravais,
                        int lmax,
                        int max_iter = 100,
                        double fstep = 0.02
                        );

  /**
   * Number of lattice vectors
   */
  int num_latt_vectors(const Matrix<double> &brav,
                       double cut,
                       const std::vector<int> &nm);


  /**
   * Get radius of truncation sphere
   */
  double trunc_radius(const Matrix<double> &brav,
                      int lmax,
                      double eta,
                      const std::vector<int> &nm);

  /**
   * Get size of lattice multiplications
   */
  std::vector<int> real_space_multiplication(const Matrix<double> &brav,
                                             int lmax,
                                             double eta);


  /**
  * Get size of reciprocal lattice multiplications
  */
  std::vector<int> reciprocal_space_multiplication(const Matrix<double> &brav,
                                             int lmax,
                                             double eta);




  double calculate_eta(Matrix<double> &brav);

  Matrix<double> calculate_dl_factor(int lmax);

  void calculate_madelung(
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
      std::vector<double> &knlatsq);


}

#endif //SRC_MADELUNG_MADELUNG_HPP
