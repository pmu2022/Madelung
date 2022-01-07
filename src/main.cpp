//
// Created by F.Moitzi on 05.01.2022.
//

#include <assert.h>

#include <any>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <variant>

#include "MultipoleMadelung.hpp"
#include "config.hpp"
#include "io.hpp"
#include "yaml-cpp/yaml.h"

int main(int argc, char *argv[]) {
  std::cout << "argc == " << argc << '\n';

  for (int ndx{}; ndx != argc; ++ndx) {
    std::cout << "argv[" << ndx << "] == " << std::quoted(argv[ndx]) << '\n';
  }
  std::cout << "argv[" << argc << "] == " << static_cast<void *>(argv[argc])
            << '\n';

  /*
   * Load a yaml script
   */
  auto file = std::string(argv[1]);

  std::cout << "File: " << file << std::endl;

  auto factory = lsms::YAMLConfigFactory(file);

  auto config = factory.createConfig();

  print_config(config);

  // input parameter
  int num_atoms = config.atom_position.n_col();

  // output parameter
  Matrix<double> madsum;
  Array3d<std::complex<double>> dl_matrix;
  Matrix<double> dl_factor;

  // Calculate the madelung matrix and all values
  std::tie(madsum, dl_matrix, dl_factor) = lsms::get_madelung(
      num_atoms, config.lmax, config.lattice, config.atom_position);

  return argc == 3 ? EXIT_SUCCESS : EXIT_FAILURE;  // optional return value
}
