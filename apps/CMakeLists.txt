add_executable(ising ising.f90 regular_graphs.f90 random_graphs.f90)
target_link_libraries(ising sigma)
