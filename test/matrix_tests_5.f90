!--------------------------------------------------------------------------!
program matrix_tests_5                                                     !
!--------------------------------------------------------------------------!
!     This program tests multiplying two sparse matrices.                  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! Matrices and graphs
    class(graph), pointer :: gr, hr, g
    integer, allocatable :: A(:,:), B(:,:), C(:,:)
    ! Graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64)
    ! Integer indices
    integer :: i,j,k,l,m,n,d,di,d1,d2,nn
    integer, allocatable :: neighbors(:), more_neighbors(:)
    ! Random numbers and vectors
    real(dp) :: p
    real(dp), allocatable :: x(:), y(:), z(:)
    ! other variables
    logical :: correct, found, isomorphic
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
    call init_seed()
    nn = 64
    p = 4.0/nn



    !----------------------------------------------------------------------!
    ! Construct reference graphs from which all test graphs are copied     !
    !----------------------------------------------------------------------!
    allocate(ll_graph::gr)
    allocate(ll_graph::hr)
    allocate(cs_graph::g)
    allocate(A(nn,nn), B(nn,nn), C(nn,nn))

    ! Make gr a ring graph
    call gr%init(nn,nn,degree=3)
    do i=1,nn
        call gr%add_edge(i,i)

        j = mod(i,nn)+1
        call gr%add_edge(i,j)
        call gr%add_edge(j,i)
    enddo

    ! Copy hr from gr
    call hr%init(gr)

    ! Compute the product of gr and hr, which should yield a ring graph
    ! where each node is connected to the next 2 nodes around the ring
    ! instead of just the next node
    call graph_product(g,gr,hr)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B)
    call hr%to_dense_graph(C)

    A = A-matmul(B,C)

    if (maxval(A)/=0) then
        print *, 'Computing graph product of ring graph with itself'
        print *, 'failed. Terminating.'
        call exit(1)
    endif

    if (verbose) print *, 'o Done testing graph product for ring graphs'


    !----------------------------------------------------------------------!
    ! Try it with different graphs                                         !
    !----------------------------------------------------------------------!
    call g%destroy()
    call gr%destroy()
    call hr%destroy()
    call gr%init(nn,nn,degree=2)

    do i=1,nn
        call gr%add_edge(i,i)

        j = mod(i,nn)+1
        call gr%add_edge(i,j)
    enddo

    call hr%init(gr,.true.)

    call graph_product(g,gr,hr)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B)
    call hr%to_dense_graph(C)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of ring graph with its reverse'
        print *, 'graph failed. Terminating.'
        call exit(1)
    endif



    !----------------------------------------------------------------------!
    ! And some more graphs                                                 !
    !----------------------------------------------------------------------!
    call g%destroy()
    call gr%destroy()
    call hr%destroy()
    call gr%init(nn,nn,degree=2)
    call hr%init(nn,nn,degree=2)

    do i=1,nn
        call gr%add_edge(i,i)
        j = mod(i,nn)+1
        call gr%add_edge(i,j)

        k = mod(i+1,nn)+1
        call hr%add_edge(i,j)
        call hr%add_edge(j,k)
    enddo

    call graph_product(g,gr,hr)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B)
    call hr%to_dense_graph(C)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of ring graph with 2-jump ring'
        print *, 'graph failed. Terminating.'
        call exit(1)
    endif

    if (verbose) print *, 'o Done testing graph product for directed graphs'



    !----------------------------------------------------------------------!
    ! And this time a random graph                                         !
    !----------------------------------------------------------------------!
    call g%destroy()
    call gr%destroy()
    call hr%destroy()
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

    call graph_product(g,gr,hr)

    if (verbose) then
        print *, 'o Done computing product g of random graphs gr, hr'
        print *, '    Number of edges:',g%ne
        print *, '    Maximum degree: ',g%max_degree
    endif

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B)
    call hr%to_dense_graph(C)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of two random graphs failed.'
        print *, 'Terminating.'
        call exit(1)
    endif


    if (verbose) print *, 'o Done testing product of random graphs'



    !----------------------------------------------------------------------!
    ! Now try it with some graph transposes                                !
    !----------------------------------------------------------------------!
    call g%destroy()

    call graph_product(g, gr, hr, trans1 = .false., trans2 = .true.)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B)
    call hr%to_dense_graph(C, trans = .true.)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of random graph with transpose'
        print *, 'of random graph failed. Terminating.'
        call exit(1)
    endif


    call g%destroy()
    call graph_product(g, gr, hr, trans1 = .true., trans2 = .false.)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B, trans = .true.)
    call hr%to_dense_graph(C)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of transpose of random graph'
        print *, 'with another random graph failed. Terminating.'
        call exit(1)
    endif


    call g%destroy()
    call graph_product(g, gr, hr, trans1 = .true., trans2 = .true.)

    call g%to_dense_graph(A)
    call gr%to_dense_graph(B, trans = .true.)
    call hr%to_dense_graph(C, trans = .true.)

    A = A-matmul(B,C)
    if (maxval(A)/=0) then
        print *, 'Computing graph product of transposes of random graphs'
        print *, 'failed. Terminating.'
        call exit(1)
    endif


    !----------------------------------------------------------------------!
    ! Clear all the graph data                                             !
    !----------------------------------------------------------------------!
    call g%destroy()
    call gr%destroy()
    call hr%destroy()


end program matrix_tests_5
