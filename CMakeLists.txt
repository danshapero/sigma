cmake_minimum_required (VERSION 2.6.2)
project( fempack Fortran C )

enable_language(Fortran)
enable_testing()

set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -fbounds-check")
set(CMAKE_Fortran_MODULE_DIRECTORY ${fempack_BINARY_DIR}/include)
set(LIBRARY_OUTPUT_PATH ${fempack_BINARY_DIR}/lib)

add_subdirectory("${CMAKE_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_SOURCE_DIR}/test")
add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
