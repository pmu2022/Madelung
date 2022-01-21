

target_sources(MultipoleMadelung PUBLIC
        Array3d.hpp
        Matrix.hpp
        NDArray.hpp
        common.hpp
        debug.hpp
        )

target_include_directories(MultipoleMadelung PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})