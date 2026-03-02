add_library(ear_star_algorithm STATIC
    input_ear_star.cpp
    ear_star_vis.cpp
    poset_ear_star.cpp
    poset_ear_star_vis.cpp
    ear_star_tri_word.cpp
)

target_include_directories(ear_star_algorithm PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}   # completion
    ${PROJECT_SOURCE_DIR}         # input.h, geometry_utils.h etc in project root
    ${PROJECT_SOURCE_DIR}/flip_algorithm # for flip algorithm
    ${PROJECT_SOURCE_DIR}/poset_visualization # for poset visualization
    
)

target_link_libraries(ear_star_algorithm PUBLIC
    polyscope
    CGAL::CGAL
)