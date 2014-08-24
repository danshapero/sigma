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
    type(ll_graph) :: g

    ! sparse and dense matrices
    class(sparse_matrix_interface), pointer :: A

    ! integer indices
    integer :: i, j, nn, frmt

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! other junk
    real(dp) :: z



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
    call g%init(nn, nn)

    do i = 1, nn
        call g%add_edge(i, i)

        j = mod(i, nn) + 1
        call g%add_edge(i, j)
        call g%add_edge(j, i)
    enddo



    !----------------------------------------------------------------------!
    ! Test each matrix type                                                !
    !----------------------------------------------------------------------!
    do frmt = 1, num_matrix_types
        if (verbose) print *, 'Format #', frmt

        ! Make a sparse matrix `A` on the reference graph `g`
        call choose_matrix_type(A, frmt)
        call A%init(nn, nn, g)

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
    enddo

    call g%destroy()

end program matrix_test_set_entry_with_realloc
