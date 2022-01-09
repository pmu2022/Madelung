

target_sources(MultipoleMadelung PUBLIC
        Array3d.hpp
        Matrix.hpp
        NDArray.hpp
        )

target_include_directories(MultipoleMadelung PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )