!--------------------------------------------------------------------------!
program matrix_tests_5                                                     !
!--------------------------------------------------------------------------!
!     This program tests multiplying two sparse matrices.                  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! Matrices and graphs
    class(graph), pointer :: gr, hr, g
    type(sparse_matrix) :: A, B, C
    real(dp), allocatable :: AD(:,:), BD(:,:), CD(:,:)
    ! Graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64)
    ! Integer indices
    integer :: i,j,k,l,m,n,d,di,d1,d2,nn
    integer, allocatable :: neighbors(:)
    ! Random numbers and vectors
    real(dp) :: p, q
    real(dp), allocatable :: x(:), y(:), z(:)
    ! other variables
    logical :: correct, found
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose


    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    select case(trim(arg))
        case("-v")
            verbose = .true.
        case("-V")
            verbose = .true.
        case("--verbose")
            verbose = .true.
    end select


    ! Initialize a random seed
    !call init_seed()
    nn = 64
    p = 4.0/nn


    !----------------------------------------------------------------------!
    ! Generate random graphs                                               !
    !----------------------------------------------------------------------!
    allocate(ll_graph::gr)
    allocate(ll_graph::hr)
    call gr%init(nn,nn)
    call hr%init(nn,nn)

    allocate(y(nn), z(nn))

    do i=1,nn
        call random_number(y)
        call random_number(z)
        do j=1,nn
            if (y(j)<p) call gr%add_edge(i,j)
            if (z(j)<p) call hr%add_edge(i,j)
        enddo
    enddo

    if (verbose) then
        print *, 'o Done generating random graphs gr, hr'
        print *, '    Number of vertices:',nn
        print *, '    Number of edges:   ',gr%ne, hr%ne
        print *, '    Maximum degree:    ',gr%max_degree, hr%max_degree
    endif


    !----------------------------------------------------------------------!
    ! Fill random matrices                                                 !
    !----------------------------------------------------------------------!
    call B%init(nn,nn,'row',gr)
    call C%init(nn,nn,'row',hr)

    allocate(neighbors(B%g%max_degree))
    do i=1,nn
        call B%g%get_neighbors(neighbors,i)

        d = B%g%degree(i)
        do k=1,d
            j = neighbors(k)

            call random_number(q)
            call B%set_value(i,j,2*q-1)
        enddo
    enddo
    deallocate(neighbors)

    allocate(neighbors(C%g%max_degree))
    do i=1,nn
        call C%g%get_neighbors(neighbors,i)

        d = C%g%degree(i)
        do k=1,d
            j = neighbors(k)

            call random_number(q)
            call C%set_value(i,j,2*q-1)
        enddo
    enddo


    !----------------------------------------------------------------------!
    ! Compute sparse matrix product                                        !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)
    call multiply_sparse_matrices(A,B,C,g)

    if (verbose) print *, 'o Done computing matrix product A = B*C'

    allocate(AD(nn,nn), BD(nn,nn), CD(nn,nn))
    call A%to_dense_matrix(AD)
    call B%to_dense_matrix(BD)
    call C%to_dense_matrix(CD)

    AD = AD-matmul(BD,CD)

    if (maxval(dabs(AD))>1.0e-14) then
        print *, 'Matrix multiplication failed.'
        print *, 'Terminating.'
        call exit(1)
    endif
    

    !----------------------------------------------------------------------!
    ! Clear all the graph data                                             !
    !----------------------------------------------------------------------!
    call g%destroy()
    call gr%destroy()
    call hr%destroy()
    deallocate(g)
    deallocate(gr)
    deallocate(hr)


end program matrix_tests_5
