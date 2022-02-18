
add_executable(testStepFunc testStepFunc.F90 readInputs.F90)

target_link_libraries(testStepFunc PUBLIC mstlib)