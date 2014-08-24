!--------------------------------------------------------------------------!
program matrix_test_sum                                                    !
!--------------------------------------------------------------------------!
! This program test explicitly adding two sparse matrices into a third     !
! matrix, as opposed to lazily forming an operator sum.                    !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs used as the matrix substrates
    type(ll_graph) :: g, h

    ! sparse and dense matrices
    class(sparse_matrix_interface), pointer :: A, B, C
    real(dp), allocatable :: AD(:,:), BD(:,:), CD(:,:)

    ! integer indices
    integer :: i, j, k, d, nn
    integer :: frmt1, frmt2, frmt3

    ! error in computing sparse matrix sum
    real(dp) :: misfit

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
    ! Make sparse matrices stored as dense arrays and form their sum       !
    !----------------------------------------------------------------------!

    nn = 48
    call init_seed()

    allocate(AD(nn, nn), BD(nn, nn), CD(nn, nn))
    AD = 0.0_dp
    BD = 0.0_dp
    CD = 0.0_dp

    call random_matrix(BD, nn, nn)
    call random_matrix(CD, nn, nn)

    AD = BD + CD



    !----------------------------------------------------------------------!
    ! Make graphs on those matrices                                        !
    !----------------------------------------------------------------------!
    call g%init(nn, nn)
    call h%init(nn, nn)

    do j = 1, nn
        do i = 1, nn
            if (BD(i, j) /= 0) call g%add_edge(i, j)
            if (CD(i, j) /= 0) call h%add_edge(i, j)
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Test each matrix type                                                !
    !----------------------------------------------------------------------!
    do frmt1 = 1, num_matrix_types
        call choose_matrix_type(B, frmt1)
        call B%init(nn, nn, g)

        do j = 1, nn
        do i = 1, nn
            z = BD(i, j)
            if (z /= 0.0_dp) call B%set_value(i, j, z)
        enddo
        enddo


        do frmt2 = 1, num_matrix_types
            call choose_matrix_type(C, frmt2)
            call C%init(nn, nn, h)

            do j = 1, nn
            do i = 1, nn
                z = CD(i, j)
                if (z /= 0.0_dp) call C%set_value(i, j, z)
            enddo
            enddo


            do frmt3 = 1, num_matrix_types
                call choose_matrix_type(A, frmt3)
                call A%set_dimensions(nn, nn)

                call sparse_matrix_sum(A, B, C)

                call A%to_dense_matrix(AD)

                misfit = maxval(dabs(BD + CD - AD))
                if (misfit > 1.0e-15) then
                    print *, 'Computing sparse matrix sum failed.'
                    print *, '||B + C - A|| =', misfit
                    print *, 'in the max-norm. Terminating.'
                    call exit(1)
                endif


                call A%destroy()
                deallocate(A)

            enddo   ! End of loop over frmt3

            call C%destroy()
            deallocate(C)

        enddo   ! End of loop over frmt2

        call B%destroy()
        deallocate(B)

    enddo   ! End of loop over frmt1


    call g%destroy()
    call h%destroy()
    



!====----------------------------------------------------------------------!
! Helper routines                                                          !
!====----------------------------------------------------------------------!
contains

!--------------------------------------------------------------------------!
subroutine random_matrix(A, m, n)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    real(dp), intent(inout) :: A(m, n)
    integer, intent(in) :: m, n
    ! local variables
    integer :: i, j, k
    real(dp) :: p, z, w

    A = 0.0_dp

    p = log(1.0_dp * m) / log(2.0_dp) / m

    do j = 1, n
        do i = 1, m
            call random_number(z)
            if (z < p) then
                call random_number(w)
                A(i, j) = w
            endif
        enddo
    enddo

    do i = 1, m
        k = count(A(i, :) /= 0)

        if (k == 0) then
            call random_number(z)
            j = int(z * n) + 1

            call random_number(w)
            A(i, j) = w
        endif
    enddo

    do j = 1, n
        k = count(A(:, j) /= 0)

        if (k == 0) then
            call random_number(z)
            i = int(z * m) + 1

            call random_number(w)
            A(i, j) = w
        endif
    enddo

end subroutine random_matrix



end program matrix_test_sum

