add_library(completion STATIC
    dedekind_poset.cpp
    dedekind_completion.cpp
    dedekind_cut_ui.cpp
    dedekind_cut.cpp
)

target_include_directories(completion PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}   # completion
    ${PROJECT_SOURCE_DIR}         # input.h, geometry_utils.h etc in project root
    ${PROJECT_SOURCE_DIR}/flip_algorithm # for flip algorithm
    ${PROJECT_SOURCE_DIR}/poset_visualization # for poset visualization

)

target_link_libraries(completion PUBLIC
    polyscope
    CGAL::CGAL
)