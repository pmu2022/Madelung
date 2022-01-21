

target_sources(MultipoleMadelung PUBLIC
        MultipoleMadelung.cpp
        MultipoleMadelung.hpp
        lattice_utils.cpp
        lattice_utils.hpp
        madelung.cpp
        madelung.hpp
        madelung_term.hpp)

target_include_directories(MultipoleMadelung PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
