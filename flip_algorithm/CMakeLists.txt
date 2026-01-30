add_library(flip_algorithm STATIC
    check_edges.cpp
    conforming.cpp
    conforming_insertion.cpp
    debug.cpp
    decomp_algo.cpp
    edge_flip_loop.cpp
    flip.cpp
    height_test.cpp
    insertion.cpp
    visualization.cpp
)

target_include_directories(flip_algorithm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}   # flip_algorithm/*.h
    ${PROJECT_SOURCE_DIR}         # input.h, geometry_utils.h etc in project root
)

target_link_libraries(flip_algorithm PUBLIC
    polyscope
    CGAL::CGAL
)

