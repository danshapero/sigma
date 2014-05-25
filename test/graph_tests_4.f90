!--------------------------------------------------------------------------!
program graph_tests_4                                                      !
!--------------------------------------------------------------------------!
!     This program is for testing whether graphs can accurately check      !
! whether or not they are symmetric, and that adding or removing edges in  !
! symmetric or asymmetric fashions preserves this property.                !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs
    class(graph), allocatable :: g, gr
    ! integer indices
    integer :: i,j,k,d,nn
    integer :: test
    integer, allocatable :: neighbors(:)
    ! random numbers
    real(dp) :: p, z
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose



    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    if (trim(arg)=="-v" .or. trim(arg)=="--verbose") then
        verbose = .true.
    endif

    ! Initialize a random seed
    nn = 128
    call init_seed()
    p = 4.0_dp/128

    ! Make a random reference graph from which graphs of all other types
    ! will be copied
    allocate(ll_graph::gr)
    call gr%init(nn)
    do i=1,nn
        do j=1,nn
            call random_number(z)
            if (i/=j .and. z<p) call gr%add_edge(i,j)
        enddo
    enddo
    allocate(neighbors(gr%max_degree))

    ! Loop through every graph type
    do test=1,4
        ! Allocate the graph
        select case(test)
            case(1)
                allocate(ll_graph::g)
                if (verbose) print *, 'Test 1, linked-list graph'
            case(2)
                allocate(coo_graph::g)
                if (verbose) print *, 'Test 2, coordinate graph'
            case(3)
                allocate(cs_graph::g)
                if (verbose) print *, 'Test 3, compressed sparse graph'
            case(4)
                allocate(ellpack_graph::g)
                if (verbose) print *, 'Test 4, ellpack graph'
        end select

        call g%init(gr)

        ! Check that the graph is properly reported as asymmetric
        call g%check_symmetry()
        if (g%symmetric) then
            print *, 'Test',test
            print *, 'Symmetry check falsely reported an asymmetric graph'
            print *, 'graph as a symmetric one, unless you happened to get'
            print *, 'lucky. Like: 2^256 lucky. Maybe run the test again to'
            print *, 'be really sure. Either way: terminating.'
            call exit(1)
        endif

        ! Symmetrize the graph
        do i=1,nn
            call gr%get_neighbors(neighbors,i)

            d = gr%degree(i)
            do k=1,d
                j = neighbors(k)
                call g%add_edge(j,i)
            enddo
        enddo

        ! Check that the graph is properly reported as symmetric now
        call g%check_symmetry()
        if (.not.g%symmetric) then
            print *, 'Test',test
            print *, 'Symmetry check falsely reported a symmetric graph'
            print *, 'as asymmetric.'
            print *, 'Terminating.'
            call exit(1)
        endif

        call g%destroy()
        deallocate(g)
    enddo


end program graph_tests_4
