!--------------------------------------------------------------------------!
program matrix_test_product                                                !
!--------------------------------------------------------------------------!
! This program tests explicitly multiplying two sparse matrices into a     !
! third matrix, as opposed to lazily forming an operator product.          !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs used as the matrix substrates
    class(graph_interface), pointer :: gr, hr, g, h

    ! sparse and dense matrices
    class(sparse_matrix_interface), pointer :: A, B, C
    real(dp), allocatable :: AD(:,:), BD(:,:), CD(:,:)

    ! indices
    integer :: i, j, k, d, nn
    integer :: frmt1, frmt2, frtm3, ordering1, ordering2, ordering3
    logical :: trans1, trans2, trans3

    ! error in computing sparse matrix product
    real(dp) :: misfit

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! other junk
    real(dp) :: z
    logical :: trans
    character(len=3) :: orientation1, orientation2, orientation3



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
    ! Make sparse matrices stored as dense arrays and form their product   !
    !----------------------------------------------------------------------!
    nn = 48
    call init_seed()

    allocate(AD(nn, nn), BD(nn, nn), CD(nn, nn))
    AD = 0.0_dp
    BD = 0.0_dp
    CD = 0.0_dp

    call random_matrix(BD, nn, nn)
    call random_matrix(CD, nn, nn)

    AD = matmul(BD, CD)



    !----------------------------------------------------------------------!
    ! Make graphs on those matrices                                        !
    !----------------------------------------------------------------------!
    allocate(ll_graph :: gr)
    allocate(ll_graph :: hr)
    call gr%init(nn, nn)
    call hr%init(nn, nn)

    do j = 1, nn
        do i = 1, nn
            if (BD(i, j) /= 0) call gr%add_edge(i, j)
            if (CD(i, j) /= 0) call hr%add_edge(i, j)
        enddo
    enddo

    call convert_graph_type(gr, "compressed sparse")
    call convert_graph_type(hr, "compressed sparse")



    !----------------------------------------------------------------------!
    ! Test each matrix type                                                !
    !----------------------------------------------------------------------!
    do frmt1 = 1, num_graph_types
    do ordering1 = 1, 2
        call choose_graph_type(g, frmt1)
        call g%copy(gr, trans)

        if (ordering1 == 1) then
            orientation1 = "row"
            trans1 = .false.
        else
            orientation1 = "col"
            trans1 = .true.
        endif

        B => sparse_matrix(nn, nn, g, orientation1)

        do j = 1, nn
        do i = 1, nn
            z = BD(i, j)
            if (z /= 0.0_dp) call B%set_value(i, j, z)
        enddo
        enddo


        do frmt2 = 1, num_graph_types
        do ordering2 = 1, 2
            call choose_graph_type(h, frmt2)
            call h%copy(hr, trans)

            if (ordering2 == 1) then
                orientation2 = "row"
                trans2 = .false.
            else
                orientation2 = "col"
                trans2 = .true.
            endif

            C => sparse_matrix(nn, nn, h, orientation2)

            do j = 1, nn
            do i = 1, nn
                z = CD(i, j)
                if (z /= 0.0_dp) call C%set_value(i, j, z)
            enddo
            enddo


            do ordering3 = 1, 2
                if (ordering3 == 1) then
                    orientation3 = "row"
                    trans3 = .false.
                else
                    orientation3 = "col"
                    trans3 = .true.
                endif

                allocate(cs_matrix :: A)
                call A%set_ordering(orientation3)
                call sparse_matrix_product(A, B, C)

                call A%to_dense_matrix(AD)

                misfit = maxval(dabs(matmul(BD, CD) - AD))
                if (misfit > 1.0e-15) then
                    print *, 'Computing sparse matrix product failed.'
                    print *, '||B * C - A|| =', misfit
                    print *, 'in the max-norm. Terminating.'
                    call exit(1)
                endif


                call A%destroy()
                deallocate(A)

            enddo

            call C%destroy()
            deallocate(C)

        enddo   ! End of loop over ordering2
        enddo   ! End of loop over frmt2

        call B%destroy()
        deallocate(B)

    enddo   ! End of loop over ordering1
    enddo   ! End of loop over frmt1




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



end program matrix_test_product

