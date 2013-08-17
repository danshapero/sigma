program solver_tests

use fempack

implicit none

    class(graph), allocatable :: g
    class(sparse_matrix), allocatable :: A
    class(iterative_solver), allocatable :: solver
    integer :: i,j,k
    real(dp) :: z
    integer, allocatable :: edges(:,:)
    real(dp), allocatable :: x(:), b(:)

    allocate(ll_graph::g)
    call g%init(99,99)
    do i=1,98
        call g%add_edge(i,i+1)
        call g%add_edge(i+1,i)
        call g%add_edge(i,i)
    enddo
    call g%add_edge(99,99)
    call convert(g,'coo')

    allocate(coo_matrix::A)
    call A%assemble(g)

    do i=1,98
        call A%set_value(i,i,2.0_dp)
        call A%set_value(i,i+1,-1.0_dp)
        call A%set_value(i+1,i,-1.0_dp)
    enddo
    call A%set_value(99,99,2.0_dp)

    allocate(x(99),b(99))
    x = 0.0_dp
    b = 2.0*(0.01_dp)**2

    allocate(cg_solver::solver)
    call solver%init(99)

    call solver%solve(A,x,b)

    print *, minval(dabs(x)), maxval(dabs(x))


end program solver_tests
