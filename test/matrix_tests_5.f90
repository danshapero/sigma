!--------------------------------------------------------------------------!
program matrix_tests_5                                                     !
!--------------------------------------------------------------------------!
!     This program tests multiplying two sparse matrices.                  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! Matrices and graphs
    class(graph), pointer :: gr, hr, g
    ! Graph edge iterators

    ! Integer indices
    integer :: i,j,k,di,nn
    integer, allocatable :: neighbors(:)
    ! Random numbers and vectors
    real(dp) :: p
    ! other variables
    logical :: correct, found
    ! command-line arguments
    character(len=16) :: arg
    logical verbose


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
    p = 6.0/nn



    !----------------------------------------------------------------------!
    ! Construct reference graphs from which all test graphs are copied     !
    !----------------------------------------------------------------------!
    allocate(ll_graph::gr)
    allocate(ll_graph::hr)
    allocate(cs_graph::g)

    call gr%init(nn,nn,degree=3)
    do i=1,nn
        call gr%add_edge(i,i)

        j = mod(i,nn)+1
        call gr%add_edge(i,j)
        call gr%add_edge(j,i)
    enddo

    call hr%init(gr)
    call graph_product(g,gr,hr)

    allocate(neighbors(g%max_degree))
    call g%get_neighbors(neighbors,16)

    correct = .true.
    do k=1,5
        found = .false.
        j = neighbors(k)

        do i=14,18
            if (j==i) found = .true.
        enddo
        correct = correct .and. found
    enddo

    if (.not.correct) then
        print *, 'Computing graph product failed, should have node 16'
        print *, 'neighboring nodes 14-18; neighbors found:',neighbors
        print *, 'Terminating.'
        call exit(1)
    endif


    !----------------------------------------------------------------------!
    ! Try it with different graphs                                         !
    !----------------------------------------------------------------------!
    call g%free()
    call gr%free()
    call hr%free()
    call gr%init(nn,nn,degree=2)

    do i=1,nn
        call gr%add_edge(i,i)

        j = mod(i,nn)+1
        call gr%add_edge(i,j)
    enddo

    call hr%init(gr,.true.)

    call graph_product(g,gr,hr)

    if (g%ne/=3*nn) then
        print *, 'Graph product does not have the correct number of edges;'
        print *, 'should be',3*nn
        print *, 'Number of edges found:',g%ne
        print *, 'Terminating.'
        call exit(1)
    endif

    if (g%max_degree/=3) then
        print *, 'Degree of graph product incorrect, should have degree = 3'
        print *, 'Degree found:',g%max_degree
        print *, 'Terminating.'
        call exit(1)
    endif

    i = 15
    do di=-1,1
        j = mod(i+di-1,nn)+1

        if (.not.g%connected(i,j)) then
            print *, 'Should have nodes',i,j
            print *, 'connected in g, but they are not!'
            print *, 'Terminating.'
            call exit(1)
        endif
    enddo


    !----------------------------------------------------------------------!
    ! And some more graphs                                                 !
    !----------------------------------------------------------------------!
    call g%free()
    call gr%free()
    call hr%free()
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



end program matrix_tests_5
