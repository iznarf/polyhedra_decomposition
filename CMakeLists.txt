cmake_minimum_required(VERSION 3.14)
project(polyhedra_decomposition)
set(CMAKE_CXX_STANDARD 17)

# Polyscope
include(FetchContent)
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG master
)
FetchContent_MakeAvailable(polyscope)

# Portable File Dialogs
FetchContent_Declare(
    portable-file-dialogs
    GIT_REPOSITORY https://github.com/samhocevar/portable-file-dialogs.git
    GIT_TAG main
)
FetchContent_MakeAvailable(portable-file-dialogs)

# vcpkg 
find_package(CGAL CONFIG REQUIRED)

add_subdirectory(flip_algorithm)

add_subdirectory(poset_visualization)

add_subdirectory(completion)

add_subdirectory(ear_star_algorithm)

add_executable(main
    main.cpp
    input.cpp
    geometry_utils.cpp
    replay.cpp
    ui_callbacks.cpp
    poset.cpp
    poset2.cpp
    compare.cpp
    interval_ui.cpp
    meet_join_ui.cpp
    moebius.cpp
    moebius_ui.cpp
    compare_nodes.cpp
    regularity_check.cpp
)


target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR})

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}/flip_algorithm
)

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}/poset_visualization
)

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}/completion
)

target_include_directories(main PRIVATE
    ${CMAKE_SOURCE_DIR}/ear_star_algorithm
)


target_include_directories(main PRIVATE ${portable-file-dialogs_SOURCE_DIR})


target_link_libraries(main PRIVATE
    flip_algorithm
    poset_visualization
    completion
    ear_star_algorithm
    polyscope
    CGAL::CGAL
)