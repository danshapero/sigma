set( test_list
type_tests.f90
util_tests.f90
graph_tests_1.f90
graph_tests_2.f90
graph_tests_3.f90
graph_tests_4.f90
matrix_test_basics.f90
matrix_test_assembly.f90
matrix_test_set_entry_with_realloc.f90
#test/matrix_tests_1.f90
#test/matrix_tests_2.f90
#test/matrix_tests_3.f90
#test/matrix_tests_4.f90
#test/matrix_tests_5.f90
#test/linear_operator_tests_1.f90
#test/solver_tests_1.f90
#test/solver_tests_2.f90
#test/wrapper_tests.c
)

#file(GLOB test_list RELATIVE ${PROJECT_SOURCE_DIR}/test ${PROJECT_SOURCE_DIR}/test/*.f90)

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
