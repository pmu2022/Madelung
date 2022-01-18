cmake_minimum_required(VERSION 3.19)
project(src/lib/Madelung)

enable_language(CXX Fortran C)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)

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
#
#
FetchContent_Declare(
        yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG        0579ae3d976091d7d664aa9d2527e0d0cff25763
)
FetchContent_MakeAvailable(yaml-cpp)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/config.yaml
        DESTINATION ${CMAKE_BINARY_DIR}/bin/)



add_subdirectory(src)

add_subdirectory(reference)


#
# Tests
#

# N. B. has to be at the top level
enable_testing()
add_subdirectory(tests)