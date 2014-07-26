!--------------------------------------------------------------------------!
program graph_test_copy                                                    !
!--------------------------------------------------------------------------!
!     This program tests whether the graph copy constructor works. A       !
! random linked-list graph `h` is generated, from which a graph `g` is     !
! copied for each graph type. The two are then checked for isomorphism.    !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs
    class(graph), pointer :: g, h

    ! dense graphs
    integer, allocatable :: A(:,:), B(:,:)

    ! integer indices
    integer :: i, j, nn, frmt, ordering

    ! random numbers
    real(dp) :: p, z

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! miscellaneous
    logical :: trans



    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!

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



    !----------------------------------------------------------------------!
    ! Set the graph size and initialize a random seed                      !
    !----------------------------------------------------------------------!

    nn = 64
    p = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Create a random sparse graph                                         !
    !----------------------------------------------------------------------!

    ! Allocate and initialize `h`
    allocate(ll_graph :: h)
    call h%init(nn, nn)

    ! Randomly add edges
    do i = 1, nn
        do j = 1, nn
            call random_number(z)

            if (z < p) call h%add_edge(i, j)
        enddo

        if (h%degree(i) == 0) then
            print *, 'Degree of node', i, 'is zero. Watch for errors in'
            print *, 'ellpack format!'
        endif
    enddo

    ! Copy `h` to a dense array `B`
    allocate(A(nn, nn), B(nn, nn))

    if (verbose) then
        print *, 'Random graph generated.'
        print *, 'Number of vertices:', nn
        print *, 'Number of edges:   ', h%ne
        print *, 'Max vertex degree: ', h%max_degree()
        print *, ' '
    endif



    !----------------------------------------------------------------------!
    ! Copy the reference graph h to graphs in other formats                !
    !----------------------------------------------------------------------!

    if (verbose) print *, 'Testing copy constructor.'

    do frmt = 1, 4
    do ordering = 1, 2
        if (verbose) print *, 'Format #', frmt, ';   ordering: ', ordering

        A = 0
        B = 0

        ! Set the variable `trans` to be true if we're testing copying the
        ! transpose of `h`
        trans = .false.
        if (ordering == 2) trans = .true.

        ! Allocate the graph `g`
        call choose_graph_type(g, frmt)

        ! Copy `g` from the reference graph `h`
        call g%copy(h, trans)

        ! Convert `g`, `h` to dense graphs
        call g%to_dense_graph(A)
        call h%to_dense_graph(B, trans)


        ! Check that the dense graphs are equal
        if ( maxval(abs(A - B)) /= 0 ) then
            print *, '    Copy constructor failed.'
            call exit(1)
        endif


        ! Deallocate `g` so it can be used for the next test
        call g%destroy()
        deallocate(g)
    enddo   ! End of loop over `ordering`
    enddo   ! End of loop over `frmt`


    deallocate(A, B)
    call h%destroy()
    deallocate(h)


end program graph_test_copy
