
target_sources(
        mstlib
        PUBLIC
        GroupCommModule.F90
        MPPModule.F90
        ParallelIOModule.F90)

target_include_directories(mstlib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR})