cmake_minimum_required(VERSION 3.19)
project(Madelung)

enable_language(CXX Fortran)

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


#
# Tests
#

enable_testing()

add_executable(legendre_test legendre.f90 legendre.h legendre_test.cpp)
target_link_libraries(legendre_test gtest_main)

# Discover all tests
include(GoogleTest)
gtest_discover_tests(legendre_test)


