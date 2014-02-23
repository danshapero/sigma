add_executable(poisson poisson.f90 fem.f90)
target_link_libraries(poisson fempack)

add_executable(permutation permutation.f90 fem.f90)
target_link_libraries(permutation fempack)

add_executable(graph_example_1 graphs/graph_example_1.f90)
target_link_libraries(graph_example_1 fempack)

add_executable(graph_example_2 graphs/graph_example_2.f90)
target_link_libraries(graph_example_2 fempack)

add_executable(graph_example_3 graphs/graph_example_3.f90)
target_link_libraries(graph_example_3 fempack)

add_executable(matrix_example_1 matrices/matrix_example_1.f90)
target_link_libraries(matrix_example_1 fempack)

add_executable(solver_example_1 solvers/solver_example_1.f90)
target_link_libraries(solver_example_1 fempack)
