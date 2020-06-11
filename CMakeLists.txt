cmake_minimum_required(VERSION 3.9)

project(lp_wrappers)
set (CMAKE_CXX_STANDARD 17)

include(cmake/FindCPLEX.cmake)
include(cmake/FindGUROBI.cmake)
find_package(Eigen3 REQUIRED)
add_subdirectory(third_party/abseil-cpp)


add_library(lp_wrappers INTERFACE)
target_include_directories(
        lp_wrappers
        INTERFACE
        include
        ${GUROBI_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIR}
)
target_link_libraries(
        lp_wrappers
        INTERFACE
        ilocplex
        ${CMAKE_DL_LIBS}
        ${GUROBI_CXX_LIBRARY}
        ${GUROBI_LIBRARY}
        absl::strings
)

add_library(
        lp_wrappers_problem
        INTERFACE
)
target_include_directories(
        lp_wrappers_problem
        INTERFACE
        include
        ${EIGEN3_INCLUDE_DIR}
)


add_executable(
    general_lp_solver
    example/general_lp_solver.cpp
)

target_link_libraries(
    general_lp_solver
    PUBLIC
    lp_wrappers
)