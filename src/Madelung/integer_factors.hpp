//
// Created by F.Moitzi on 24.12.2021.
//

#ifndef MADELUNG_INTEGER_FACTORS_HPP
#define MADELUNG_INTEGER_FACTORS_HPP

#include <vector>

namespace lsms {

  std::vector<int> get_lofk(int l0);

  std::vector<int> get_mofk(int l0);

  std::vector<int> get_lofj(int l0);

  std::vector<int> get_mofj(int l0);

  std::vector<int> get_kofj(int l0);

}

#endif //MADELUNG_INTEGER_FACTORS_HPP
