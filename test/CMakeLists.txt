set( test_list
type_tests.f90
util_tests.f90
graph_test_basics.f90
graph_test_copy.f90
matrix_test_basics.f90
matrix_test_set_multiple_entries.f90
matrix_test_set_entry_with_realloc.f90
matrix_test_sum.f90
matrix_test_product.f90
matrix_test_strategy.f90
linear_operator_test_algebra.f90
solver_test_diffusion_1d.f90
solver_test_advection_diffusion_1d.f90
solver_test_jacobi.f90
solver_test_incomplete_cholesky.f90
#test/matrix_tests_4.f90
#test/matrix_tests_5.f90
#test/wrapper_tests.c
)


foreach(test ${test_list})
    string(REPLACE ".f90" "" test ${test})
    add_executable(${test} ${PROJECT_SOURCE_DIR}/test/${test}.f90)
    target_link_libraries(${test} sigma)
    add_test(${test} ${test})
endforeach(test)

#set(C_HEADER_FILES ${PROJECT_SOURCE_DIR}/include/graphs.h)
#add_executable(wrapper_tests ${PROJECT_SOURCE_DIR}/test/wrapper_tests.c ${C_HEADER_FILES})
#target_link_libraries(wrapper_tests sigma)
#add_test(wrapper_tests wrapper_tests)
