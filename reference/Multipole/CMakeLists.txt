


#
# Testing purposes
#

add_library(multipole_madelung SHARED)

target_sources(
        multipole_madelung
        PUBLIC
        multipole_madelung.h
        IntegerFactorsModule.f90
        KindParamModule.f90
        MadelungModule.f90
        ErrorHandlerModule.f90
        GauntFactorsModule.f90
        LegendreModule.f90
        SphericalHarmonicsModule.f90
        MathParamModule.f90
        LatticeModule.f90
        gaussq.f90
        breakLine.f90
        CMadelungModule.f90
        ord3v.f90
)

target_compile_options(
        multipole_madelung
        PUBLIC
        "-cpp"
)

if (NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR "Generic")
endif ()
find_package(LAPACK REQUIRED)
target_link_libraries(multipole_madelung PUBLIC LAPACK::LAPACK)

find_package(BLAS REQUIRED)
target_link_libraries(multipole_madelung PUBLIC BLAS::BLAS)

target_include_directories(
        multipole_madelung
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)