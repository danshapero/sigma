!--------------------------------------------------------------------------!
program matrix_test_speed_query                                            !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph representing the matrix structure
    class(graph_interface), pointer :: g

    ! sparse matrices
    type(sparse_matrix) :: L
    class(sparse_matrix_interface), pointer :: A

    ! integer indices
    integer :: i, j, k, nn

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
    ! Set the matrix size and initialize a random seed                     !
    !----------------------------------------------------------------------!
    nn = 256
    probability = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random graph                                                  !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        do j = i + 1, nn
            call random_number(z)
            if (z < probability) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Test that basic sparse matrix formats recognize whether getting a    !
    ! row or column is fast                                                !
    !----------------------------------------------------------------------!

    allocate(csr_matrix::A)
    call graph_laplacian(A, g)

    if (.not. A%is_get_row_fast()) call fail("row", "CSR", "fast")
    if (A%is_get_column_fast()) call fail("column", "CSR", "slow")

    call A%destroy()
    deallocate(A)


    allocate(csc_matrix::A)
    call graph_laplacian(A, g)

    if (A%is_get_row_fast()) call fail("row", "CSC", "slow")
    if (.not. A%is_get_column_fast()) call fail("col", "CSC", "fast")

    call A%destroy()
    deallocate(A)


    allocate(ellpack_matrix::A)
    call graph_laplacian(A, g)

    if (.not. A%is_get_row_fast()) call fail("row", "ellpack", "fast")
    if (A%is_get_column_fast())    call fail("col", "ellpack", "slow")

    call A%destroy()
    deallocate(A)


    allocate(default_row_matrix::A)
    call graph_laplacian(A, g)

    if (.not. A%is_get_row_fast()) call fail("row", "default row", "fast")
    if (A%is_get_column_fast())    call fail("col", "default row", "slow")

    call A%destroy()
    deallocate(A)


    allocate(default_column_matrix::A)
    call graph_laplacian(A, g)

    if (A%is_get_row_fast()) call fail("row", "default column", "slow")
    if (.not. A%is_get_column_fast()) call fail("col", "default column", &
                                                                    & "fast")

    call A%destroy()
    deallocate(A)



    !----------------------------------------------------------------------!
    ! Test that the composite sparse matrix recognizes whether getting a   !
    ! row or column is fast when used as a wrapper around a basic type.    !
    !----------------------------------------------------------------------!

    call L%set_matrix_type("csr")

    call graph_laplacian(L, g)

    if (.not. L%is_get_row_fast()) then
        call fail("row", "wrapper around CSR", "fast")
        call exit(1)
    endif

    if (L%is_get_column_fast()) then
        call fail("column", "wrapper around CSR", "slow")
        call exit(1)
    endif

    call L%destroy()
    call g%destroy()



    !----------------------------------------------------------------------!
    ! Test that the composite sparse matrix recognizes whether getting a   !
    ! row or column is fast when used as a composite of several matrices.  !
    !----------------------------------------------------------------------!

    call L%set_dimensions(nn, nn)
    call L%set_block_sizes([nn/2, nn/2], [nn/2, nn/2])

    call g%init(nn/2, nn/2)
    probability = 2 * log(0.5_dp * nn) / log(2.0_dp) / nn
    do i = 1, nn/2
        call g%add_edge(i, i)
        do j = i + 1, nn/2
            call random_number(z)
            if (z < probability) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    allocate(csr_matrix::A)
    call A%init(nn/2, nn/2, g)
    call L%set_submatrix(1, 1, A)

    call A%remove_reference()
    nullify(A)

    allocate(ellpack_matrix::A)
    call A%init(nn/2, nn/2, g)
    call L%set_submatrix(2, 2, A)

    call A%remove_reference()
    nullify(A)

    allocate(csr_matrix::A)
    call A%init(nn/2, nn/2, g)
    call L%set_submatrix(1, 2, A)

    call A%remove_reference()
    nullify(A)

    allocate(default_row_matrix::A)
    call A%init(nn/2, nn/2, g)
    call L%set_submatrix(2, 1, A)

    call A%remove_reference()
    nullify(A)


    if (.not. L%is_get_row_fast()) call fail("row", "composite row", "fast")
    if (L%is_get_column_fast()) call fail("col", "composite row", "slow")


    allocate(csc_matrix::A)
    call A%init(nn/2, nn/2, g)
    call L%set_submatrix(2, 1, A)

    call A%remove_reference()
    nullify(A)


    if (L%is_get_row_fast()) call fail("row", "mixed composite", "slow")
    if (L%is_get_column_fast()) call fail("col", "mixed composite", "slow")


    call L%destroy()
    call g%destroy()
    deallocate(g)

contains


!--------------------------------------------------------------------------!
subroutine graph_laplacian(A, g)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(graph_interface), intent(in) :: g
    ! local variables
    integer :: i, j, k, d, nn
    integer, allocatable :: nodes(:)

    d = g%get_max_degree()
    allocate(nodes(d))

    nn = g%n

    call A%init(nn, nn, g)

    do i = 1, nn
        call g%get_neighbors(nodes, i)
        d = g%get_degree(i)

        do k = 1, d
            j = nodes(d)

            call A%add(i, j, -1.0_dp)
            call A%add(i, i, +1.0_dp)
        enddo
    enddo

end subroutine graph_laplacian



!--------------------------------------------------------------------------!
subroutine fail(row_or_col, matrix_type, fast_or_not)                      !
!--------------------------------------------------------------------------!
    character(len=*), intent(in) :: row_or_col, matrix_type, fast_or_not

    print *, "Getting ", row_or_col, " of ", matrix_type, &
        & " matrix should be ", fast_or_not
    call exit(1)

end subroutine fail



end program matrix_test_speed_query


