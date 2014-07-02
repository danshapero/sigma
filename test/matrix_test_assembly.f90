!--------------------------------------------------------------------------!
program matrix_test_assembly                                               !
!--------------------------------------------------------------------------!
! This program tests the `assemble` method for matrices, which renders the !
! underlying graph structure immutable and makes for faster matrix-vector  !
! multiplication.                                                          !
!--------------------------------------------------------------------------!

    use sigma

    implicit none

    ! graph used as the matrix substrate
    class(graph), pointer :: g, h

    ! sparse and dense matrices
    class(sparse_matrix), pointer :: A
    real(dp), allocatable :: B(:,:), C(:,:)

    ! vectors
    real(dp), allocatable :: x(:), y1(:), y2(:)

    ! integer indices
    integer :: i, j, k, d, nn, frmt, ordering, ind(2), ord(2)

    ! random numbers
    real(dp) :: p, w, z

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! other junk
    character(len=3) :: orientation
    integer :: row_degree, col_degree



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
    ! Set the matrix size and initialize a random seed                     !
    !----------------------------------------------------------------------!

    nn = 64
    p = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random reference sparse matrix, stored as a dense matrix      !
    !----------------------------------------------------------------------!

    allocate(B(nn, nn))
    allocate(x(nn), y1(nn), y2(nn))
    B = 0.0_dp

    do i = 1, nn
        do j = 1, nn
            call random_number(z)
            if (z < p) then
                call random_number(w)
                B(i, j) = w
            endif
        enddo
    enddo

    ! Compute the column degree of the graph
    row_degree = 0
    do i = 1, nn
        k = count(B(i, :) /= 0)
        row_degree = max(k, row_degree)
    enddo

    ! Compute the row degree of the graph
    col_degree = 0
    do j = 1, nn
        k = count(B(:, j) /= 0)
        col_degree = max(k, col_degree)
    enddo

    d = max(row_degree, col_degree)



    !----------------------------------------------------------------------!
    ! Test each matrix type                                                !
    !----------------------------------------------------------------------!
    do frmt = 1, 4
    do ordering = 1, 2
        orientation = "row"
        ord = [1, 2]
        if (ordering == 2) then
            orientation = "col"
            ord = [2, 1]
        endif

        if (verbose) print *, 'Format #',frmt, '; order: ',orientation

        call choose_graph_type(g, frmt)

        ! Initialize `g` so that it will have a fixed amount of extra space
        ! that we can compress away later
        call g%init(nn, nn, degree = d + 3)

        ! ...which means we have to copy from the dense matrix manually
        do j = 1, nn
            do i = 1, nn
                if (B(i, j) /= 0) then
                    ind = [i, j]
                    ind = ind(ord)
                    call g%add_edge(ind(1), ind(2))
                endif
            enddo
        enddo

        ! Initialize `A` with the unnecessarily large graph
        A => sparse_matrix(nn, nn, g, orientation)

        ! Set the entries of `A` to be the same as those of `B`
        do j = 1, nn
            do i = 1, nn
                z = B(i, j)
                if (z /= 0) call A%set_value(i, j, z)
            enddo
        enddo

        !--------
        ! Assemble `A`
        call A%assemble()


        !--------
        ! Check that the entries of `A` are still the same as those of `B`
        do j = 1, nn
            do i = 1, nn
                z = B(i, j)
                if (A%get_value(i, j) /= z) then
                    print *, i, j, z, A%get_value(i, j)
                    print *, 'Assembly of A failed, matrix entries have '
                    print *, 'changed. Terminating.'
                    call exit(1)
                endif
            enddo
        enddo


        !--------
        ! Check that matrix multiplication works
        call random_number(x)
        y1 = 0.0_dp
        y2 = 0.0_dp

        call A%matvec(x, y1)
        y2 = matmul(B, x)

        w = maxval(dabs(y1 - y2)) / maxval(dabs(y2))
        if (w > 1.0e-15) then
            print *, 'Matrix-vector multiplication gives wrong result'
            print *, 'after matrix assembly. Terminating.'
            call exit(1)
        endif

    enddo
    enddo



end program matrix_test_assembly
