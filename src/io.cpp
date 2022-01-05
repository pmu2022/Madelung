//
// Created by F.Moitzi on 05.01.2022.
//

#include "io.hpp"

#include "yaml-cpp/yaml.h"

#include "Matrix.hpp"


lsms::MadelungConfig lsms::YAMLConfigFactory::createConfig() {

  YAML::Node yaml = YAML::LoadFile("config.yaml");

  // L cutoff
  auto lmax = yaml["lmax"].as<int>();

  // Name
  auto name = yaml["name"].as<std::string>();

  // Lattice
  Matrix<double> lattice(3, 3);
  YAML::Node node = yaml["lattice"];
  for (std::size_t i = 0; i < node.size(); i++) {
    YAML::Node subnode = node[i];
    for (std::size_t j = 0; j < subnode.size(); j++) {
      lattice(i, j) = subnode[j].as<double>();
    }
  }

  // Atom positions
  YAML::Node atom_node = yaml["atomposition"];

  auto natoms = atom_node.size();
  Matrix<double> atom_position(3, natoms);

  for (std::size_t i = 0; i < natoms; i++) {
    for (std::size_t j = 0; j < 3; j++) {
      atom_position(j, i) = atom_node[i][j].as<double>();
    }
  }


  return {name, lmax, lattice, atom_position};
}

lsms::YAMLConfigFactory::YAMLConfigFactory(const std::string &path) : path{path} {

}
