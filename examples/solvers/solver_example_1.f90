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
    class(graph), pointer :: g
    real(dp) :: p
    real(dp), allocatable :: z(:)

    ! a sparse matrix
    class(sparse_matrix), pointer :: A

    ! a solver object
    type(cg_solver) :: solver
    type(jacobi_solver) :: pc

    ! some integer indices
    integer :: i, j, k, nn

    ! variables for making the graph Laplacian of g
    integer :: d
    integer, allocatable :: neighbors(:)

    ! some vectors
    real(dp), allocatable :: x(:), b(:)



    nn = 1024

    ! initialize a random seed
    call init_seed()
    p = 5.0 / nn

    allocate(x(nn), b(nn), z(nn))



    !----------------------------------------------------------------------!
    ! Set up a random graph g                                              !
    !----------------------------------------------------------------------!
    ! Make h a linked-list graph, since it's cheap to modify this type
    allocate(ll_graph::g)
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        call random_number(z)

        do j = i + 1, nn
            if (z(j) < p) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    ! Convert `g` to a nicer format
    call convert_graph_type(g, 'compressed sparse')



    !----------------------------------------------------------------------!
    ! Make A the graph Laplacian of g + the identity matrix                !
    !----------------------------------------------------------------------!
    A => sparse_matrix(nn, nn, g, 'row')
    call A%zero()

    d = g%max_degree()
    allocate(neighbors(d))

    ! For each vertex of the graph,
    do i = 1, nn
        ! set A[i,i] = 1.0.
        call A%set_value(i, i, 1.0_dp)

        d = g%degree(i)
        call g%get_neighbors(neighbors, i)

        ! For each neighbor j of i,
        do k = 1, d
            j = neighbors(k)

            if (j /= i) then
                ! Set A[i,j] = -1.0, and A[i,i] = A[i,i]+1.0.
                call A%set_value(i, j, -1.0_dp)
                call A%add_value(i, i, 1.0_dp)
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
    call solver%solve(A, x, b, pc)

    ! Report some results.
    write(*,100) solver%iterations
100 format('Conjugate gradient solver completed in ', i4, ' steps.')
    write(*,120) minval(x), maxval(x)
120 format('Range of solution: (', e12.6, ', ', e12.6, ').')



    call g%destroy()
    call A%destroy()
    deallocate(g)
    deallocate(x,b)


end program solver_example_1
