//
// Created by F.Moitzi on 16.12.2021.
//

#ifndef MADELUNG_MULTIPOLEMADELUNG_HPP
#define MADELUNG_MULTIPOLEMADELUNG_HPP

#include "Matrix.hpp"
#include "Array3d.hpp"

#include <vector>
#include <complex>


namespace lsms {

  constexpr auto EPSI = 1e-14;



  std::tuple<Matrix<double>, Array3d<std::complex<double>>, Matrix<double> > calculate_reduced(int num_atoms,
                                                                                            int lmax,
                                                                                            Matrix<double> &bravais,
                                                                                            Matrix<double> &atom_position);


  double calculate_eta(Matrix<double> &brav);

  Matrix<double> calculate_dl_factor(int kmax_mad, int jmax_mad, int lmax_mad);

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


  void init_madelung(Matrix<double> &bravais,
                     Matrix<double> &atom_position,
                     int lmax);


  double scaling_factor(Matrix<double> &bravais, int lmax);

  void reciprocal_lattice(Matrix<double> &bravais,
                          Matrix<double> &reciprocal_bravais,
                          double &scale);

  void real_space_trunc(Matrix<double> &brav,
                        int lmax_mad,
                        double eta,
                        double &rscut,
                        std::vector<int> &nm);

  void reciprocal_space_trunc(Matrix<double> &brav,
                              int lmax_mad,
                              double eta,
                              double &rscut,
                              std::vector<int> &nm);


  __attribute__((unused)) int num_latt_vectors_legacy(Matrix<double> &brav,
                                                      double cut,
                                                      std::vector<int> &nm);

  int num_latt_vectors(Matrix<double> &brav,
                       double cut,
                       std::vector<int> &nm);


  void generate_lattice(Matrix<double> &r_brav,
                        Matrix<double> &k_brav,
                        int lmax,
                        double eta);

}


#endif //MADELUNG_MULTIPOLEMADELUNG_HPP
