add_executable(poisson poisson.f90 fem.f90)
target_link_libraries(poisson fempack)

add_executable(permutation permutation.f90 fem.f90)
target_link_libraries(permutation fempack)

add_executable(graph_example_1 graphs/graph_example_1.f90)
target_link_libraries(graph_example_1 fempack)
