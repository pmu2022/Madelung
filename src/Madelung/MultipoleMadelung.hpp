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



  class MultipoleMadelung {

  private:

    /// Global number of atoms
    int num_atoms;

    /// Local number of atoms
    int local_num_atoms;

    /// Madelung sum
    Matrix<double> madsum;

  public:

    /// Lattice of the structure
    const Matrix<double> &lattice;

    /// Atomic positions
    Matrix<double> atom_position;

    /// Angular-momentum index cutoff l
    int lmax;

    /// Index of local atoms
    const std::vector<int> &global_position_index;


    Array3d<std::complex<double>> dl_matrix;

    Matrix<double> dl_factor;

    MultipoleMadelung(const Matrix<double> &lattice,
                      Matrix<double> atom_position,
                      int lmax,
                      const std::vector<int> &global_position_index);

    double getMadSum(int i, int j) const;



  };



}


#endif //MADELUNG_MULTIPOLEMADELUNG_HPP
