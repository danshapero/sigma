cmake_minimum_required (VERSION 2.6.2)
project( sigma Fortran C )

enable_language(Fortran)
enable_language(CXX)
enable_testing()

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -fbounds-check")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    set (CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -check -fpe0")
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${sigma_BINARY_DIR}/include)
set(LIBRARY_OUTPUT_PATH ${sigma_BINARY_DIR}/lib)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
endif()

include_directories("${PROJECT_SOURCE_DIR}/include")

add_subdirectory("${CMAKE_SOURCE_DIR}/src")
add_subdirectory("${CMAKE_SOURCE_DIR}/test")
add_subdirectory("${CMAKE_SOURCE_DIR}/examples")
add_subdirectory("${CMAKE_SOURCE_DIR}/apps")
