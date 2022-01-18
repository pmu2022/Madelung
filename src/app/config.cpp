//
// Created by F.Moitzi on 05.01.2022.
//

#include "config.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <variant>

#include "common.hpp"

template <class T>
static std::ostream &operator<<(std::ostream &os, const lsms::matrix<T> &v) {
  for (auto i = 0; i < v.shape()[0]; i++) {
    for (auto j = 0; j < v.shape()[1]; j++) {
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
  std::cout << "Atom position: " << std::endl
            << config.atom_position << std::endl;
}
