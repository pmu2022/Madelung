
#add_library(IntegrationVP SHARED)

find_package(MPI REQUIRED)

add_library(mstlib SHARED)

add_subdirectory(lib)

add_subdirectory(plib)

add_subdirectory(src)

target_include_directories(mstlib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(mstlib
        PUBLIC
        MPI::MPI_Fortran
        )

target_compile_definitions(mstlib
        PUBLIC
        MPI
        )


target_compile_options(mstlib
        PUBLIC
        "-fallow-argument-mismatch")

add_subdirectory(driver)

#add_subdirectory(src)

