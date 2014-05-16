!--------------------------------------------------------------------------!
program graph_tests_3                                                      !
!--------------------------------------------------------------------------!
!     This program performs further tests for adding new edges to a graph, !
! in the exceptional case where a graph no longer has any storage space    !
! and re-allocation of the graph's internal structure may be necessary.    !
!--------------------------------------------------------------------------!

use sigma

implicit none

    class(graph), allocatable :: g
    integer :: i,j,k,d,nn,test
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose


    ! Get command line arguments to see if we're running in verbose mode
    verbose = .false.
    call getarg(1,arg)
    select case(trim(arg))
        case("-v","-V","--verbose")
            verbose = .true.
    end select

    nn = 128

    do test=1,4
        ! Allocate and initialize the graph
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

        call g%init(nn,degree=2)

        ! Make g a ring graph, where each vertex is connected to the next
        ! and previous vertex mod 128.
        do i=1,nn
            j = mod(i,nn)+1
            call g%add_edge(i,j)

            j = mod(i+nn-2,nn)+1
            call g%add_edge(i,j)
        enddo

        if (verbose) print *, '    Done creating 2-regular ring graph.'

        do i=1,nn
            if (g%degree(i)/=2) then
                print *, 'Test',test
                print *, 'Should have a regular graph where degree of'
                print *, 'every vertex is',2
                print *, 'Degree of',i,'is',g%degree(i)
                call exit(1)
            endif
        enddo


        ! Add a bunch more edges that the graph doesn't have room for,
        ! and see if the program fails
        d = int(log(1.0_dp*nn)/log(2.0_dp))
        do i=1,nn
            do k=2,d
                j = mod(i+k,nn)+1
                call g%add_edge(i,j)

                j = mod(i+127-k,nn)+1
                call g%add_edge(i,j)
            enddo
        enddo

        if (verbose) print *, '    Done adding more edges to ring graph.'

        ! Check if the degrees of the nodes have been updated properly
        do i=1,nn
            if (g%degree(i)/=2*d) then
                print *, 'Test',test
                print *, 'Should have a regular graph where degree of'
                print *, 'every vertex is',2*d
                print *, 'Degree of',i,'is',g%degree(i)
                print *, 'Terminating.'
                call exit(1)
            endif
        enddo

        ! Check if the edges are connected properly
        do i=1,nn
            do k=2,d
                j = mod(i+k,nn)+1
                if (.not.g%connected(i,j)) then
                    print *, 'Test',test
                    print *, 'Edges',i,j,'not connected, even though they'
                    print *, 'were added to g! Terminating.'
                    call exit(1)
                endif

                j = mod(i+127-k,nn)+1
                if (.not.g%connected(i,j)) then
                    print *, 'Test',test
                    print *, 'Edges',i,j,'not connected, even though they'
                    print *, 'were added to g! Terminating.'
                    call exit(1)
                endif
            enddo
        enddo

        call g%free()
        deallocate(g)
        if (verbose) print *, ' '
    enddo

end program graph_tests_3
