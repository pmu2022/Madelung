//
// Created by F.Moitzi on 15.12.2021.
//

#ifndef SPHERICAL_HARMONICS_HPP
#define SPHERICAL_HARMONICS_HPP

#include <complex>

#include <vector>
#include <array>

namespace lsms {

  namespace math {

    std::vector<std::complex<double>> spherical_harmonics(std::array<double, 3> vec,
                                                          unsigned lmax);

  }

}


#endif //SPHERICAL_HARMONICS_HPP
