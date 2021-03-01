cmake_minimum_required(VERSION 3.9)

project(lp_wrappers)
set (CMAKE_CXX_STANDARD 17)

option(LPWRAPPERS_WITH_CPLEX "with cplex" ON)
option(LPWRAPPERS_WITH_GUROBI "with gurobi" ON)
option(LPWRAPPERS_BUILD_EXAMPLES "build examples" ON)

if(LPWRAPPERS_WITH_CPLEX)
    include(cmake/FindCPLEX.cmake)
    SET(LPWRAPPERS_CPLEX_LIBRARIES ilocplex)
endif()

if(LPWRAPPERS_WITH_GUROBI)
    include(cmake/FindGUROBI.cmake)
endif()

find_package(Eigen3 REQUIRED)

if(NOT TARGET absl::strings)
    add_subdirectory(third_party/abseil-cpp)
endif()

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
        ${LPWRAPPERS_CPLEX_LIBRARIES}
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

if(LPWRAPPERS_BUILD_EXAMPLES)
    if(LPWRAPPERS_WITH_CPLEX)
        add_executable(
            cplex_lp_multiple_times
            example/cplex_lp_multiple_times.cpp
        )

        target_link_libraries(
            cplex_lp_multiple_times
            PUBLIC
            lp_wrappers
        )
    endif()
    if (LPWRAPPERS_WITH_GUROBI)
        add_executable(
            gurobi_lp_multiple_times
            example/gurobi_lp_multiple_times.cpp
        )

        target_link_libraries(
            gurobi_lp_multiple_times
            PUBLIC
            lp_wrappers
        )
    endif()
endif()