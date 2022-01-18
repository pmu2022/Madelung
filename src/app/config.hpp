//
// Created by F.Moitzi on 05.01.2022.
//

#ifndef MADELUNG_CONFIG_HPP
#define MADELUNG_CONFIG_HPP

#include <any>
#include <map>
#include <string>
#include <utility>

#include "common.hpp"

namespace lsms {

class MadelungConfig {
 public:
  MadelungConfig(const std::string &name,
                 int lmax,
                 const lsms::matrix<double> &lattice,
                 const lsms::matrix<double> &atom_position)
      : name{name},
        lmax{lmax},
        lattice{lattice},
        atom_position{atom_position} {};

  std::string name;
  int lmax;
  matrix<double> lattice;
  matrix<double> atom_position;
};

void print_config(MadelungConfig &config);

}  // namespace lsms

#endif  // MADELUNG_CONFIG_HPP
