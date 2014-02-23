program solver_tests

use fempack

implicit none

    class(graph), pointer :: g
    type(sparse_matrix) :: A
    class(iterative_solver), allocatable :: solver
    class(preconditioner), allocatable :: pc
    integer :: i,n,test,nnz_per_row(99)
    real(dp) :: c,dx
    real(dp), allocatable :: u(:), b(:), u_c(:)

    nnz_per_row = 3
    allocate(coo_graph::g)
    call g%init(99,99,nnz_per_row)
    do i=1,98
        call g%add_edge(i,i+1)
        call g%add_edge(i+1,i)
        call g%add_edge(i,i)
    enddo
    call g%add_edge(99,99)

    call A%init(99,99,'row',g)

    do i=1,98
        call A%set_value(i,i,2.0_dp)
        call A%set_value(i,i+1,-1.0_dp)
        call A%set_value(i+1,i,-1.0_dp)
    enddo
    call A%set_value(99,99,2.0_dp)

    allocate(u(99),b(99),u_c(99))
    dx = 0.01_dp
    b = 2.0*dx**2

    do n=1,99
        u_c(n) = n*dx*(1.0_dp-n*dx)
    enddo

    do test=1,2
        u = 0.0_dp

        select case(test)
            case(1)
                allocate(cg_solver::solver)
                print *, 'CG solver test'
            case(2)
                allocate(bicgstab_solver::solver)
                print *, 'BiCG solver test'
        end select

        call solver%init(99,tolerance=1.0d-16)

        allocate(jacobi_preconditioner::pc)
        call pc%init(A,0)

        call solver%solve(A,u,b,pc)

        if ( maxval(dabs(u-u_c))>1.0e-14 ) then
            print *, 'Max value of u should be 1/4;'
            print *, 'Value found:', maxval(u)
            call exit(1)
        endif

        deallocate(solver,pc)
    enddo

    ! Tests for non-symmetric matrices
    c = 0.5_dp
    do n=1,98
        call A%set_value(n,n,2.0_dp)
        call A%set_value(n+1,n,-(1.0_dp+c*dx))
        call A%set_value(n,n+1,-(1.0_dp-c*dx))
    enddo
    call A%set_value(99,99,2.0_dp)

    do n=1,99
        u_c(n) = exp(0.5*c*n*dx)*sin(pi*n*dx)
    enddo
    call A%matvec(u_c,b)

    allocate(bicgstab_solver::solver)
    allocate(jacobi_preconditioner::pc)
    call solver%init(99,tolerance=1.0d-14)
    call pc%init(A,0)

    u = 0.0_dp

    call solver%solve(A,u,b,pc)

    if( maxval(dabs(u-u_c))>1.0e-12 ) then
        print *, 'BiCG-Stab solver failed for non-symmetric matrix'
        print *, maxval(dabs(u-u_c)),maxval(dabs(u))
        call exit(1)
    endif



end program solver_tests
