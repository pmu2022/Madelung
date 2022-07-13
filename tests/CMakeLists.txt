

#
# Tests
#

# Copy test files

add_custom_target(copy-test-files ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/RhoTildaL)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/RhoTildaL
        COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/data/RhoTildaL
        ${CMAKE_CURRENT_BINARY_DIR}/RhoTildaL
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/RhoTildaL)



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

add_executable(
        regression_supercell_test regression_supercell_test.cpp
)
target_link_libraries(regression_supercell_test PUBLIC multipole_madelung)
target_link_libraries(regression_supercell_test PUBLIC MultipoleMadelung)
target_link_libraries(regression_supercell_test PUBLIC gtest_main)

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
target_link_libraries(boost_test PUBLIC Boost::boost)

# 8. Math libary wrapper
add_executable(
        math_test math_test.cpp
)
target_link_libraries(math_test PUBLIC multipole_madelung)
target_link_libraries(math_test PUBLIC MultipoleMadelung)
target_link_libraries(math_test PUBLIC gtest_main)

# 9. Multipole moments tests
add_executable(
        multipole_moments_tests multipole_moments_tests.cpp
)
target_link_libraries(multipole_moments_tests PUBLIC multipole_madelung)
target_link_libraries(multipole_moments_tests PUBLIC MultipoleMadelung)
target_link_libraries(multipole_moments_tests PUBLIC gtest_main)
add_dependencies(multipole_moments_tests copy-test-files)

# Discover all tests
include(GoogleTest)
gtest_discover_tests(legendre_test)
gtest_discover_tests(spherical_harmonics_test)
gtest_discover_tests(madelung_test)
gtest_discover_tests(gaunt_factor_test)
gtest_discover_tests(regression_test)
gtest_discover_tests(regression_supercell_test)
gtest_discover_tests(multi_array_test)
gtest_discover_tests(boost_test)
gtest_discover_tests(math_test)
