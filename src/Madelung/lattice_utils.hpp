//
// Created by F.Moitzi on 21.12.2021.
//

#ifndef MADELUNG_LATTICE_HPP
#define MADELUNG_LATTICE_HPP

#include "Matrix.hpp"
#include "Array3d.hpp"

#include <vector>
#include <tuple>
#include <complex>

namespace lsms {

  /**
   * Create all vectors in certain cutoff with a certain repeotition
   */
  Matrix<double> create_lattice(Matrix<double> &brav,
                                double cutoff,
                                const std::vector<int> &nm,
                                int size);

  /**
    * Create all vectors in certain cutoff with a certain repeotition
    */
  std::tuple<Matrix<double>, std::vector<double>> create_lattice_and_sq(Matrix<double> &brav,
                                                                        double cutoff,
                                                                        const std::vector<int> &nm,
                                                                        int size);

  /**
   *  inserts a vector in a list of vectors such that they are in
   *  order of increasing length.
   */
  void insert_ordered(Matrix<double> &latt_vec,
                      std::vector<double> &latt_vec_sq,
                      int len,
                      std::vector<double> &vec,
                      double &v_sq);

  /**
   * Lattice volumes
   */
   double omega(Matrix<double> &bravais);

   /**
    * Calculate reciprocal lattice
    */
   void reciprocal_lattice(Matrix<double> &bravais,
                           Matrix<double> &reciprocal_bravais,
                           double &scale);

}


#endif //MADELUNG_LATTICE_HPP
