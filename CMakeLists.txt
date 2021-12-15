cmake_minimum_required(VERSION 3.19)
project(Madelung)

enable_language(CXX Fortran C)

set(CMAKE_CXX_STANDARD 14)

#
# Google tests
#
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_subdirectory(Madelung)

add_subdirectory(lib)

add_subdirectory(include)


#
# Tests
#

enable_testing()


#
# Tests all routines in the Madelung library
#

# 1. Test associated legendre polynoms
add_executable(legendre_test legendre_test.cpp)
target_link_libraries(legendre_test PUBLIC multipole_madelung)
target_link_libraries(legendre_test PUBLIC MultipoleMadelung)
target_link_libraries(legendre_test PUBLIC gtest_main)

# 2. Test spherical harmonics
add_executable(spherical_harmonics_test
        spherical_harmonics_test.cpp
)
target_link_libraries(spherical_harmonics_test PUBLIC multipole_madelung)
target_link_libraries(spherical_harmonics_test PUBLIC MultipoleMadelung)
target_link_libraries(spherical_harmonics_test PUBLIC gtest_main)

# 3. Test Madelung potentials
add_executable(
        madelung_test madelung_test.cpp
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

# Discover all tests
include(GoogleTest)
gtest_discover_tests(legendre_test)
gtest_discover_tests(spherical_harmonics_test)
gtest_discover_tests(gaunt_factor_test)
gtest_discover_tests(madelung_test)



