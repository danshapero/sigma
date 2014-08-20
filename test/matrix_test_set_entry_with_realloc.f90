!--------------------------------------------------------------------------!
program matrix_test_set_entry_with_realloc                                 !
!--------------------------------------------------------------------------!
! This program tests setting a matrix entry (i, j) which is not already an !
! edge of the underlying graph. This operation necessitates reallocating   !
! some of the matrix structure, which should be avoided.                   !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    class(graph_interface), pointer :: g, h

    ! sparse and dense matrices
    class(sparse_matrix_interface), pointer :: A

    ! integer indices
    integer :: i, j, nn, frmt, ordering

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! other junk
    real(dp) :: z
    logical :: trans
    character(len=3) :: orientation



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
    ! Make a reference graph                                               !
    !----------------------------------------------------------------------!

    nn = 64

    allocate(ll_graph :: h)
    call h%init(nn, nn)

    do i = 1, nn
        call h%add_edge(i, i)

        j = mod(i, nn) + 1
        call h%add_edge(i, j)
        call h%add_edge(j, i)
    enddo



    !----------------------------------------------------------------------!
    ! Test each matrix type                                                !
    !----------------------------------------------------------------------!
    do frmt = 1, 4
    do ordering = 1, 2
        orientation = "row"
        trans = .false.

        if (ordering == 2) then
            orientation = "col"
            trans = .true.
        endif

        if (verbose) print *, 'Format #',frmt,'; order: ',orientation

        ! Make a copy `g` of `h`, possibly with edges reversed
        call choose_graph_type(g, frmt)
        call g%copy(h, trans)

        ! Make a sparse matrix `A` on that graph
        A => sparse_matrix(nn, nn, g, orientation)

        ! Set the entries of `A`
        do i = 1, nn
            call A%set_value(i, i, +2.0_dp)

            j = mod(i, nn) + 1
            call A%set_value(i, j, -1.0_dp)
            call A%set_value(j, i, -1.0_dp)
        enddo

        !--------
        ! Set entries of `A` that have not been pre-allocated in the graph
        do i = 1, nn
            j = mod(i + 1, nn) + 1
            call A%set_value(i, j, -1.0_dp)
            call A%add_value(i, i, +1.0_dp)
        enddo

        !--------
        ! Check that the entries have been set properly
        do i = 1, nn
            j = mod(i + 1, nn) + 1
            z = A%get_value(i, j)

            if (z /= -1) then
                print *, 'Matrix entry not set. Terminating.'
                call exit(1)
            endif
        enddo


        call A%destroy()
        deallocate(A)
        call g%destroy()
        deallocate(g)
    enddo
    enddo



end program matrix_test_set_entry_with_realloc
