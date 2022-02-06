//
// Created by F.Moitzi on 17.12.2021.
//

#ifndef MADELUNG_UTILS_HPP
#define MADELUNG_UTILS_HPP

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>

template<typename Iter_T>
auto norm(Iter_T first, Iter_T last) {
  return std::sqrt(std::inner_product(first, last, first, 0.0));
}

template<typename Iter_T>
auto norm_sq(Iter_T first, Iter_T last) {
  return std::inner_product(first, last, first, 0.0);
}

namespace lsms {

  namespace math {

    template<typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
    inline constexpr T pow(const T base, unsigned const exponent) {
      return (exponent == 0) ? 1 : (base * pow(base, exponent - 1));
    }

    template<typename T>
    std::vector<double> linspace(T start_in, T end_in, std::size_t num_in) {

      std::vector<double> linspaced;

      double start = static_cast<double>(start_in);
      double end = static_cast<double>(end_in);
      double num = static_cast<double>(num_in);

      if (num == 0) { return linspaced; }
      if (num == 1) {
        linspaced.push_back(start);
        return linspaced;
      }

      double delta = (end - start) / (num - 1);

      for (int i = 0; i < num - 1; ++i) {
        linspaced.push_back(start + delta * i);
      }

      // Ensure that start and end are exactly the same as the input
      linspaced.push_back(end);
      return linspaced;
    }

  }

}

/**
 *
 *  call calGammaFunc to calculate the integral:
 *
 *         inf       2*l
 *  I(l) = int dx * x   * exp(-x**2)
 *          a
 *
 *  and store it in gamma_l(l)
 *
 *  l = 0, 1, ..., lmax.
 *
 *  using: (2*l+1)*I(l) = 2*I(l+1) - c(l) ,
 *
 *
 *         c(l) = a**(2*l+1) * exp(-a**2)
 *
 *         I(0) = sqrt(pi)/2 * erfc(a)
 *
 *  erfc(z) = 1 - erf(z)
 */
template<class T>
std::vector<T> gamma_func(const T &&a, int &lmax) {

  auto r_factor = 0.5 * std::exp(-a * a);

  std::vector<T> gamma(lmax + 1, 0.0);

  // l = 0
  gamma[0] = std::sqrt(M_PI) * 0.5 * std::erfc(a);

  for (int l = 1; l <= lmax; l++) {
    gamma[l] = 0.5 * (2 * l - 1) * gamma[l - 1] + std::pow(a, (2 * l - 1)) * r_factor;
  }

  return gamma;

}

#endif //MADELUNG_UTILS_HPP
