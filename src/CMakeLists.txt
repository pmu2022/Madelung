
add_library(MultipoleMadelung SHARED)

add_subdirectory(common)

add_subdirectory(Madelung)


add_executable(Madelung main.cpp)

set_target_properties(Madelung PROPERTIES OUTPUT_NAME madelung)

set_target_properties(Madelung
        PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        )

target_sources(Madelung PUBLIC
        config.cpp
        config.hpp
        io.cpp
        io.hpp
        )

target_link_libraries(Madelung PUBLIC MultipoleMadelung)


target_link_libraries(Madelung PUBLIC yaml-cpp)
target_include_directories(Madelung PUBLIC ${yaml-cpp_SOURCE_DIR}/include)