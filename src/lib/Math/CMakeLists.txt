



target_sources(MultipoleMadelung PUBLIC
        clm.f90
        CMakeLists.txt
        gaunt_factor.cpp
        gaunt_factor.f90
        gaunt_factor.hpp
        gaussq.f90
        gaussq.h
        integer_factors.cpp
        integer_factors.f90
        integer_factors.hpp
        legendre.f90
        legendre.h
        legendre.hpp
        integrator.hpp
        utils.hpp
        spherical_harmonics.f90
        spherical_harmonics.hpp
        spherical_harmonics.cpp
        spherical_harmonics.h
        utils.hpp
        legendre.hpp)


target_include_directories(MultipoleMadelung PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})