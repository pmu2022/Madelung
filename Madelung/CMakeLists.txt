

add_library(MultipoleMadelung SHARED)

target_sources(MultipoleMadelung
        PUBLIC
        MultipoleMadelung.cpp
        MultipoleMadelung.hpp
        Lattice.cpp
        Lattice.hpp
        Matrix.hpp
        Array3d.hpp
        utils.hpp
        debug.hpp
        common.hpp
        spherical_harmonics.f90
        spherical_harmonics.hpp
        integer_factors.f90
        legendre.f90
        clm.f90
        integer_factors.hpp
        gaunt_factor.f90
        gaunt_factor.hpp
        gaussq.f90
        gaunt_factor.cpp
        )

target_include_directories(MultipoleMadelung PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )