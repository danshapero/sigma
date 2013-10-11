add_executable(poisson poisson.f90 fem.f90)
target_link_libraries(poisson fempack)

add_executable(permutation permutation.f90 fem.f90)
target_link_libraries(permutation fempack)
