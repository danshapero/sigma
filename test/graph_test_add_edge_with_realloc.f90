!--------------------------------------------------------------------------!
program graph_test_add_edge_with_realloc                                   !
!--------------------------------------------------------------------------!
!     This program tests adding new edges to a graph, in the exceptional   !
! case where a graph no longer has any storage space and reallocation of   !
! the graph's internal structure may be necessary.                         !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph
    class(graph), pointer :: g

    ! integer indices
    integer :: i, j, k, d, nn, frmt

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose



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


    nn = 128

    do frmt = 1, 4
        if (verbose) print *, 'Format #', frmt

        ! Allocate and initialize the graph
        call choose_graph_type(g, frmt)
        call g%init(nn, degree = 2)


        ! Make `g` a ring graph, where each vertex is connected to the next
        ! and previous vertex around the ring modulo `nn`
        do i = 1, nn
            j = mod(i, nn) + 1
            call g%add_edge(i, j)
            call g%add_edge(j, i)
        enddo

        if (verbose) print *, '    Done creating 2-regular ring graph.'


        ! Check that it was set right
        do i = 1, nn
            if (g%degree(i) /= 2) then
                print *, 'Should have a regular graph where degree of'
                print *, 'every vertex is ', 2
                print *, 'Degree of vertex', i, 'is', g%degree(i)
            endif
        enddo


        ! Add more edges, for which the graph doesn't have room.
        d = int( log(1.0_dp * nn) / log(2.0_dp) )

        do i = 1, nn
            do k = 2, d
                j = mod(i + k, nn) + 1
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            enddo
        enddo

        if (verbose) print *, '    Done adding more edges to ring graph.'


        ! Check that the number of edges has been updated properly
        if (g%ne /= 2 * d * nn) then
            print *, 'After adding edges to graph, total number of edges'
            print *, 'should be:', 2 * d * nn
            print *, 'Number of edges found:', g%ne
        endif


        ! Check that the correct degree of each vertex is returned
        do i = 1, nn
            if (g%degree(i) /= 2 * d) then
                print *, 'Should have a regular graph where degree of'
                print *, 'every vertex is',2*d
                print *, 'Degree of',i,'is',g%degree(i)
                print *, 'Terminating.'
                call exit(1)
            endif
        enddo


        ! Check that the connectivity structure is right
        do i = 1, nn
            do k = 2, d
                j = mod(i + k, nn) + 1
                if (.not. (g%connected(i, j) .and. g%connected(j, i))) then
                    print *, 'Edges', i, j, 'not connected, even though they'
                    print *, 'were added to g! Terminating.'
                    call exit(1)
                endif
            enddo
        enddo


        ! Clear the graph data structure for the next test
        call g%destroy()
        deallocate(g)
        if (verbose) print *, ' '
    enddo   ! End of loop over `frmt`



end program graph_test_add_edge_with_realloc
