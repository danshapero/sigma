!--------------------------------------------------------------------------!
program matrix_test_copy                                                   !
!--------------------------------------------------------------------------!
!     This program tests copying one matrix into another. A random sparse  !
! matrix `B` is generated in each format, which is then copied into a      !
! sparse matrix `A` of each format. The two are checked for equality in    !
! every element.                                                           !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs
    type(ll_graph) :: g

    ! sparse matrices
    class(sparse_matrix_interface), pointer :: A, B, C

    ! integer indices
    integer :: i, j, k, d, nn, frmt1, frmt2
    integer, allocatable :: nodes(:)

    ! random numbers
    real(dp) :: probability, z

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



    !----------------------------------------------------------------------!
    ! Set the graph size and initialize a random seed                      !
    !----------------------------------------------------------------------!

    nn = 128
    probability = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Create a random sparse graph                                         !
    !----------------------------------------------------------------------!

    call g%init(nn, nn / 2)

    do i = 1, nn
        do j = 1, nn / 2
            call random_number(z)
            if (z < probability) call g%add_edge(i, j)
        enddo
    enddo

    d = g%get_max_degree()
    allocate(nodes(d))

    print *, g%get_num_edges()
    print *, " "



    !----------------------------------------------------------------------!
    ! Copy graphs to / from each type and check for validity               !
    !----------------------------------------------------------------------!

    do frmt1 = 1, num_matrix_types
        if (verbose) print *, frmt1

        call choose_matrix_type(A, frmt1)
        call A%init(nn, nn / 2, g)

        do i = 1, nn
            d = g%get_degree(i)
            call g%get_neighbors(nodes, i)

            do k = 1, d
                j = nodes(k)
                call random_number(z)
                call A%set(i, j, z)
            enddo
        enddo

        do frmt2 = 1, num_matrix_types
            if (verbose) print *, frmt2

            !-----------------------------------
            ! Check that copying a matrix works
            call choose_matrix_type(B, frmt2)
            call B%set_dimensions(nn, nn / 2)

            call B%copy_matrix(A, trans = .false.)

            do i = 1, nn
                do j = 1, nn / 2
                    if ( dabs(A%get(i, j) - B%get(i, j)) > 1.0e-14 ) then
                        print *, "Copying matrix failed", i, j
                        call exit(1)
                    endif
                enddo
            enddo

            call B%destroy()
            deallocate(B)


            !----------------------------------------------------
            ! Check that copying the transpose of a matrix works
            call choose_matrix_type(C, frmt2)
            call C%set_dimensions(nn / 2, nn)

            call C%copy_matrix(A, trans = .true.)

            do i = 1, nn
                do j = 1, nn / 2
                    if ( dabs(A%get(i, j) - C%get(j, i)) > 1.0e-14 ) then
                        print *, "Copying transpose matrix failed", i, j
                        call exit(1)
                    endif
                enddo
            enddo

            call C%destroy()
            deallocate(C)
        enddo

        call A%destroy()
        deallocate(A)
    enddo


    call g%destroy()


end program matrix_test_copy
