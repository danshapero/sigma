add_executable(ising ising.f90 regular_graphs.f90 random_graphs.f90)
target_link_libraries(ising sigma)

add_executable(self_avoiding_walk self_avoiding_walk.f90 regular_graphs.f90 random_graphs.f90)
target_link_libraries(self_avoiding_walk sigma)
