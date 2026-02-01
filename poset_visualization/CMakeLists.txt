add_library(poset_visualization STATIC
    poset_vis.cpp
    poset_vis_ui.cpp
    render_helpers.cpp
    poset_tri_mesh.cpp
    poset_vis_helpers.cpp
)

target_include_directories(poset_visualization PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}   # poset visualization 
    ${PROJECT_SOURCE_DIR}         # input.h, geometry_utils.h etc in project root
    ${PROJECT_SOURCE_DIR}/flip_algorithm # for flip algorithm
    ${PROJECT_SOURCE_DIR}/completion # for completion
)

target_link_libraries(poset_visualization PUBLIC
    polyscope
    CGAL::CGAL
)