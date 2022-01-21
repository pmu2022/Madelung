
#include "spherical_harmonics.hpp"

#include "integer_factors.hpp"

#ifdef USE_BOOST

  #include <boost/math/special_functions/spherical_harmonic.hpp>

#elif

  #include "spherical_harmonics.h"

#endif

std::vector<std::complex<double>> lsms::math::spherical_harmonics(std::array<double, 3> vec, unsigned int lmax) {


  std::vector<std::complex<double>> ylm(lsms::get_kmax(lmax));

#ifdef USE_BOOST


  auto radius = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  auto theta = std::acos(vec[2] / radius);
  auto phi = std::atan2(vec[1], vec[0]);

  auto k = 0;

  for (auto l = 0; l <= lmax; l++) {

    for (auto m = -l; m <= l; m++) {

      ylm[k] = boost::math::spherical_harmonic(l, m, theta, phi);
      k++;


    }

  }

#elif

  sph_harm_1(vec.data(), &lmax_mad, ylm.data());

#endif

  return ylm;

}
