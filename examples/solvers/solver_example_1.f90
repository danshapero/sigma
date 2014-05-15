!==========================================================================!
program solver_example_1                                                   !
!==========================================================================!
!==== Example program demonstrating the use of an iterative solver     ====!
!==== without a preconditioner.                                        ====!
!==== It is assumed that the user has some familiarity with matrix     ====!
!==== operations in order to understand this program.                  ====!
!==== As in previous examples, the graph is generated randomly.        ====!
!==========================================================================!

use sigma

implicit none

    ! a linked-list graph and variables for generating it randomly
    class(graph), pointer :: g, h
    real(dp) :: p
    real(dp), allocatable :: z(:)

    ! a sparse matrix
    type(sparse_matrix) :: A

    ! a solver object
    type(cg_solver) :: solver
    type(jacobi_solver) :: pc

    ! some integer indices
    integer :: i,j,k,nn

    ! variables for making the graph Laplacian of g
    integer :: d
    integer, allocatable :: neighbors(:)

    ! some vectors
    real(dp), allocatable :: x(:), b(:)



    nn = 1024

    ! initialize a random seed
    call init_seed()
    p = 5.0/nn

    allocate(x(nn),b(nn),z(nn))



    !----------------------------------------------------------------------!
    ! Set up a random graph g                                              !
    !----------------------------------------------------------------------!
    ! Make h a linked-list graph, since it's cheap to modify this type
    allocate(ll_graph::h)
    call h%init(nn)

    do i=1,nn
        call h%add_edge(i,i)

        call random_number(z)

        do j=i+1,nn
            if (z(j)<p) then
                call h%add_edge(i,j)
                call h%add_edge(j,i)
            endif
        enddo
    enddo

    ! Make g a compressed sparse graph, since these have faster accesses
    allocate(cs_graph::g)

    ! Copy the connectivity structure of h to g
    call g%copy(h)

    ! Compress g
    call g%compress()

    ! Delete h
    deallocate(h)



    !----------------------------------------------------------------------!
    ! Make A the graph Laplacian of g + the identity matrix                !
    !----------------------------------------------------------------------!
    call A%init(nn,nn,'row',g)
    call A%zero()

    allocate(neighbors(g%max_degree))

    ! For each vertex of the graph,
    do i=1,nn
        ! set A[i,i] = 1.0.
        call A%set_value(i,i,1.0_dp)

        d = g%degree(i)
        call g%get_neighbors(neighbors,i)

        ! For each neighbor j of i,
        do k=1,d
            j = neighbors(k)

            if (j/=0 .and. j/=i) then
                ! Set A[i,j] = -1.0, and A[i,i] = A[i,i]+1.0.
                call A%set_value(i,j,-1.0_dp)
                call A%add_value(i,i,1.0_dp)
            endif
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Make a random right-hand side, set up the iterative solver object    !
    ! and solve the system.                                                !
    !----------------------------------------------------------------------!
    call random_number(b)

    ! Initialize the solver for a system with nn unknowns and a solution
    ! tolerance of 1.0e-10.
    call solver%setup(A)
    call pc%setup(A)

    ! Call the solver on A, with x as the approximate solution and b as 
    ! the right-hand side.
    call solver%solve(A,x,b,pc)

    ! Report some results.
    write(*,100) solver%iterations
100 format('Conjugate gradient solver completed in ',i4,' steps.')
    write(*,120) minval(x),maxval(x)
120 format('Range of solution: (',e12.6,', ',e12.6,').')


end program solver_example_1
