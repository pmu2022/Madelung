

#
# Tests
#

find_package(GSL REQUIRED)

# 1. Test associated legendre polynoms
add_executable(legendre_test legendre_test.cpp)
target_link_libraries(legendre_test PUBLIC multipole_madelung)
target_link_libraries(legendre_test PUBLIC MultipoleMadelung)
target_link_libraries(legendre_test PUBLIC gtest_main)
target_link_libraries(legendre_test PUBLIC GSL::gsl)

# 2. Test spherical harmonics
add_executable(spherical_harmonics_test
        spherical_harmonics_test.cpp
        )
target_link_libraries(spherical_harmonics_test PUBLIC multipole_madelung)
target_link_libraries(spherical_harmonics_test PUBLIC MultipoleMadelung)
target_link_libraries(spherical_harmonics_test PUBLIC gtest_main)
target_link_libraries(spherical_harmonics_test PUBLIC GSL::gsl)

# 3. Test Gaunt coeffiencts
add_executable(madelung_test
        madelung_test.cpp
        )
target_link_libraries(madelung_test PUBLIC multipole_madelung)
target_link_libraries(madelung_test PUBLIC MultipoleMadelung)
target_link_libraries(madelung_test PUBLIC gtest_main)

# 4. Test Gaunt coeffiencts
add_executable(gaunt_factor_test
        gaunt_factor_test.cpp
        )
target_link_libraries(gaunt_factor_test PUBLIC multipole_madelung)
target_link_libraries(gaunt_factor_test PUBLIC MultipoleMadelung)
target_link_libraries(gaunt_factor_test PUBLIC gtest_main)

# 5. Test Madelung potentials
add_executable(
        regression_test regression_test.cpp
)
target_link_libraries(regression_test PUBLIC multipole_madelung)
target_link_libraries(regression_test PUBLIC MultipoleMadelung)
target_link_libraries(regression_test PUBLIC gtest_main)

# 6. Test Madelung potentials
add_executable(
        multi_array_test multi_array_test.cpp
)
target_link_libraries(multi_array_test PUBLIC multipole_madelung)
target_link_libraries(multi_array_test PUBLIC MultipoleMadelung)
target_link_libraries(multi_array_test PUBLIC gtest_main)

# 7. Boost
add_executable(
        boost_test boost_test.cpp
)
target_link_libraries(boost_test PUBLIC multipole_madelung)
target_link_libraries(boost_test PUBLIC MultipoleMadelung)
target_link_libraries(boost_test PUBLIC gtest_main)
target_link_libraries(boost_test PUBLIC Boost::math)
target_link_libraries(boost_test PUBLIC Boost::boost)

# Discover all tests
include(GoogleTest)
gtest_discover_tests(legendre_test)
gtest_discover_tests(spherical_harmonics_test)
gtest_discover_tests(madelung_test)
gtest_discover_tests(gaunt_factor_test)
gtest_discover_tests(regression_test)
gtest_discover_tests(multi_array_test)
gtest_discover_tests(boost_test)
