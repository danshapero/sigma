!--------------------------------------------------------------------------!
program matrix_test_basics                                                 !
!--------------------------------------------------------------------------!
! This program performs tests of basic matrix operations:                  !
!     o initialization                                                     !
!     o getting / setting entries                                          !
!     o iteration over entries                                             !
!     o matrix-vector multiplication                                       !
!     o storage compression                                                !
!     o permutation                                                        !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    class(graph), pointer :: g, h

    ! sparse and dense matrices
    class(sparse_matrix), pointer :: A
    real(dp), allocatable :: B(:,:)

    ! vectors
    real(dp), allocatable :: x(:), y(:)

    ! integer indices
    integer :: i, j, k, d, nn, frmt, ordering

    ! permutation
    integer, allocatable :: p(:)

    ! variables for getting matrix rows / columns
    integer :: row_degree, col_degree
    integer, allocatable :: nodes(:)
    real(dp), allocatable :: slice(:)

    ! random numbers
    real(dp) :: c, w, z

    ! other junk
    logical :: correct, trans
    character(len=3) :: orientation


    nn = 64
    c = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random reference sparse matrix, stored as a dense matrix      !
    !----------------------------------------------------------------------!

    allocate(B(nn, nn))
    B = 0.0_dp

    do i = 1, nn
        do j = 1, nn
            call random_number(z)
            if (z < c) then
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
    allocate( nodes(d), slice(d) )



    !----------------------------------------------------------------------!
    ! Make a graph on the reference matrix                                 !
    !----------------------------------------------------------------------!

    allocate(ll_graph :: h)
    call h%init(nn, nn)

    do i = 1, nn
        do j = 1, nn
            if (B(i,j) /= 0) call h%add_edge(i, j)
        enddo
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

            ! Make a copy `g` of `h`, possibly with the edges reversed
            call choose_graph_type(g, frmt)
            call g%copy(h, trans)

            ! Make a sparse matrix `A` on that graph
            A => sparse_matrix(nn, nn, g, orientation)

            ! Check that all the entries of `A` are zero
            do i = 1, nn
                do j = 1, nn
                    z = A%get_value(i, j)
                    if (z /= 0) then
                        print *, 'Entries of a just-initialized sparse '
                        print *, 'matrix should be zero! Terminating.'
                        call exit(1)
                    endif
                enddo
            enddo

            ! Set the entries of `A` to be the same as those of `B`
            do j = 1, nn
                do i = 1, nn
                    z = B(i, j)
                    if (z /= 0) call A%set_value(i, j, z)
                enddo
            enddo

            ! Check that the entries of `A` are the same as those of `B`
            do j = 1, nn
                do i = 1, nn
                    z = B(i, j)
                    if (A%get_value(i, j) /= z) then
                        print *, 'Setting entry of A failed. Terminating.'
                        call exit(1)
                    endif
                enddo
            enddo

            ! Check that getting an entire row / column of the matrix works
            do i = 1, nn
                call A%get_row(nodes, slice, i)

                ! First, check that every entry returned from the sparse
                ! matrix `A` corresponds to the right value in the
                ! reference matrix `B`
                do k = 1, d
                    j = nodes(k)
                    if (j /= 0) then
                        if (slice(k) /= B(i, j)) then
                            print *, 'Getting row of sparse matrix failed.'
                            print *, 'Terminating.'
                            call exit(1)
                        endif
                    endif
                enddo

                ! Next, check that every non-zero entry in row `i` of
                ! `B` was actually returned by `get_row`
                do j = 1, nn
                    if (B(i, j) /= 0) then
                        correct = .false.
                        do k = 1, d
                            correct = correct .or. (nodes(k) == j)
                        enddo

                        if (.not. correct) then
                            print *, 'Getting row of sparse matrix failed,'
                            print *, 'did not return entry that is in '
                            print *, 'the row. Terminating.'
                            call exit(1)
                        endif
                    endif
                enddo
            enddo

            do j = 1, nn
                call A%get_column(nodes, slice, j)

                do k = 1, d
                    i = nodes(k)
                    if (i /= 0) then
                        if (slice(k) /= B(i, j)) then
                            print *, 'Getting col of sparse matrix failed.'
                            print *, 'Terminating.'
                            call exit(1)
                        endif
                    endif
                enddo

                do i = 1, nn
                    if (B(i, j) /= 0) then
                        correct = .false.
                        do k = 1, d
                            correct = correct .or. (nodes(k) == i)
                        enddo

                        if (.not. correct) then
                            print *, 'Getting col of sparse matrix failed,'
                            print *, 'did not return entry that is in '
                            print *, 'the col. Terminating.'
                            call exit(1)
                        endif
                    endif
                enddo
            enddo


            ! Destroy the matrix and graph so they're ready for the next
            ! test
            call A%destroy()
            call g%destroy()
            deallocate(A)
            deallocate(g)
        enddo
    enddo


end program matrix_test_basics
