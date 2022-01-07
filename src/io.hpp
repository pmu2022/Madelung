//
// Created by F.Moitzi on 05.01.2022.
//

#ifndef MADELUNG_IO_HPP
#define MADELUNG_IO_HPP

#include <string>

#include "config.hpp"

namespace lsms {

class ConfigFactory {
 public:
  virtual MadelungConfig createConfig() = 0;
};

class YAMLConfigFactory : public ConfigFactory {
 private:
  std::string path;

 public:
  explicit YAMLConfigFactory(const std::string &path);

  MadelungConfig createConfig() override;
};

}  // namespace lsms

#endif  // MADELUNG_IO_HPP
