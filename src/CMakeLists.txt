
# Shared libaray
add_library(MultipoleMadelung SHARED)

if(USE_BOOST)
    target_compile_definitions(MultipoleMadelung PUBLIC USE_BOOST)
    target_link_libraries(MultipoleMadelung PUBLIC Boost::boost)
endif()

add_subdirectory(lib)

add_subdirectory(app)


