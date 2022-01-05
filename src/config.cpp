//
// Created by F.Moitzi on 05.01.2022.
//

#include "config.hpp"

#include "Matrix.hpp"

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include <variant>


template<class T>
static std::ostream &operator<<(std::ostream &os, Matrix<T> &v) {

  for (auto i = 0; i < v.n_row(); i++) {
    for (auto j = 0; j < v.n_col(); j++) {
      os << v(i, j) << " ";
    }
    os << std::endl;
  }

  return os;
}


void lsms::print_config(lsms::MadelungConfig &config) {

  std::cout << "Name:          " << config.name << std::endl;
  std::cout << "Lmax:          " << config.lmax << std::endl;
  std::cout << "Lattice:       " << std::endl << config.lattice << std::endl;
  std::cout << "Atom position: " << std::endl << config.atom_position << std::endl;

}
