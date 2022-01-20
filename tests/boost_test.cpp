


#include <gtest/gtest.h>

#include <complex>

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

using namespace boost::multiprecision;
using namespace boost::math;

cpp_int boost_factorial(int num) {
  cpp_int fact = 1;
  for (int i = num; i > 1; --i)
    fact *= i;
  return fact;
}

TEST(BoostTest, SphericalHarmonics) {

  unsigned int l = 1;
  int m = -1;

  std::complex<double> result = boost::math::spherical_harmonic(l, m, 0.1, 0.2);
}

TEST(BoostTest, SphericalHarmonics) {

  int num = 30;
  std::cout << "Factorial of " << num << " = " << boost_factorial(num) << std::endl;


}

