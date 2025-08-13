cmake_minimum_required(VERSION 3.10)
project(2D_DIRAC_MG)

find_package(MFEM REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Spectra)

option(BARYCENTRIC_DUAL "Use barycenters instead of circumcenters to assemble dual grid" OFF)
option(COARSE_SOLVE_PINV "Use a sparse QR factorization to compute the least-squareas solution instead of sparse LU" ON)

if(BARYCENTRIC_DUAL)
  add_definitions(-DBARYCENTRIC_DUAL)
endif()

if(COARSE_SOLVE_PINV)
  add_definitions(-DCOARSE_SOLVE_PINV)
endif()

set(CMAKE_CXX_FLAGS "-Ofast -march=native -Wall -fopenmp")

include_directories(inc)
add_subdirectory(src)
add_subdirectory(targets)
