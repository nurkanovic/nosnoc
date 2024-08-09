cmake_minimum_required(VERSION 3.10.2)
project(nosnoc_solver-nlp-sol)

# TODO this should be generalized
list(APPEND CMAKE_MODULE_PATH ~/tools/casadi/cmake)
find_package(CASADI REQUIRED)
include_directories(~/tools/casadi/build/install/include)
find_package(IPOPT REQUIRED)
link_directories(${IPOPT_LIBRARY_DIRS})

add_library(nosnoc_solver_nlp SHARED nosnoc_solver_nlp.c)

set_target_properties(nosnoc_solver_nlp PROPERTIES PREFIX "")
set_target_properties(nosnoc_solver_nlp PROPERTIES SUFFIX ".casadi")

add_library(nosnoc_solver_comp SHARED nosnoc_solver_comp.c)

set_target_properties(nosnoc_solver_comp PROPERTIES PREFIX "")
set_target_properties(nosnoc_solver_comp PROPERTIES SUFFIX ".casadi")

set(CASADI_DEPENDENCIES ${CMAKE_DL_LIBS})

message("${CASADI_LIBARARIES}")
add_executable(homotopy_solver main.cpp HomotopySolver.cpp)
target_link_libraries(homotopy_solver ${CASADI_LIBRARIES} ${CASADI_DEPENDENCIES} ${IPOPT_LIBRARIES})
target_include_directories(homotopy_solver PRIVATE .)