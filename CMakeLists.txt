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

add_executable(main
    main.cpp
    input.cpp
    visualization.cpp
    check_edges.cpp
    conforming.cpp
    flip.cpp
    insertion.cpp
    edge_flip_loop.cpp
    conforming_insertion.cpp
    debug.cpp
    geometry_utils.cpp
    replay.cpp
    ui_callbacks.cpp
    poset.cpp
)

target_include_directories(main PRIVATE ${portable-file-dialogs_SOURCE_DIR})
target_link_libraries(main PRIVATE polyscope CGAL::CGAL)