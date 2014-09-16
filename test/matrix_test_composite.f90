!--------------------------------------------------------------------------!
program matrix_test_composite                                              !
!--------------------------------------------------------------------------!
!     This program tests using the type sparse_matrix as a composite of    !
! several other sparse matrices, i.e. a block matrix. Each of the sub-     !
! matrices of the composite could be one of the primitive matrix formats,  !
! like CSC or ellpack, or could itself be a composite matrix.              !
!     Using sparse_matrix objects in this fashion allows the user to       !
! cleanly generate independent sub-blocks of a large matrix (possibly in   !
! parallel) and lace them together later.                                  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    class(graph_interface), pointer :: g

    ! sparse matrix objects
    type(sparse_matrix) :: A, B
    class(sparse_matrix_interface), pointer :: C

    ! vectors
    real(dp), allocatable :: x(:), y(:), z(:)

    ! integer indices
    integer :: i, j, k, d, nn, nn1, nn2

    ! variables for getting matrix rows / columns
    integer, allocatable :: nodes(:)

    ! random numbers
    real(dp) :: p, q

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
    ! Set the problem dimensions                                           !
    !----------------------------------------------------------------------!
    nn1 = 768
    nn2 = 512

    nn = nn1 + nn2



    !----------------------------------------------------------------------!
    ! Initialize a sparse matrix composite                                 !
    !----------------------------------------------------------------------!
    call A%set_dimensions(nn, nn)
    call A%set_block_sizes( [nn1, nn2], [nn1, nn2] )

    if (verbose) then
        print *, "o Composite sparse matrix initialized."
    endif



    !----------------------------------------------------------------------!
    ! Make a random Laplacian matrix                                       !
    !----------------------------------------------------------------------!

    ! Make the graph for the (1, 1)-submatrix of `A`
    p = log(1.0_dp * nn1) / log(2.0_dp) / nn1

    call choose_graph_type(g, "ll")

    call erdos_renyi_graph(g, nn1, nn1, p, symmetric = .true.)
    d = g%max_degree()

    if (verbose) then
        print *, "o Done generating first Erdos-Renyi graph."
        print *, "    Number of vertices: ", nn1
        print *, "    Number of edges:    ", g%ne
        print *, "    Maximum degree:     ", d
    endif

    ! Make a Laplacian matrix on the graph with random weights
    call choose_matrix_type(C, "csr")
    call erdos_renyi_matrix(C, g)

    if (verbose) then
        print *, "o Done generating first random weighted Laplacian matrix."
    endif

    ! Set the (1, 1)-submatrix of `A` to be `C`
    call A%set_submatrix(1, 1, C)

    ! Nullify `C` so that we can use it to build another sub-matrix
    call C%remove_reference()
    nullify(C)
    call g%destroy()


    ! Do the same thing again for the (2, 2)-submatrix of `A`
    p = log(1.0_dp * nn2) / log(2.0_dp) / nn2

    call erdos_renyi_graph(g, nn2, nn2, p, symmetric = .true.)
    d = g%max_degree()

    if (verbose) then
        print *, "o Done generating second Erdos-Renyi graph."
        print *, "    Number of vertices: ", nn2
        print *, "    Number of edges:    ", g%ne
        print *, "    Maximum degree:     ", d
    endif

    call choose_matrix_type(C, "csr")
    call erdos_renyi_matrix(C, g)

    if (verbose) then
        print *, "o Done generating second random weighted Laplacian matrix."
    endif

    call A%set_submatrix(2, 2, C)

    call C%remove_reference()
    nullify(C)
    call g%destroy()


    ! Now generate a graph for the (1, 2)-submatrix
    p = 6.0 / nn1
    call erdos_renyi_graph(g, nn1, nn2, p, symmetric = .false.)

    ! Convert it to CS storage
    call convert_graph_type(g, "cs")

    ! Make the (1, 2)-submatrix of `A` a CSR matrix,
    call A%set_matrix_type(1, 2, "csr")

    ! then *set* its graph to point to `g`.
    call A%set_graph_submat(1, 2, g)

    ! Now make the (2, 1)-submatrix of `A` a CSC matrix,
    call A%set_matrix_type(2, 1, "csc")

    ! and set its graph to point to `g` also.
    call A%set_graph_submat(2, 1, g)

    ! There are now 3 references to `g`: one from us having created it in
    ! the first place, another from the (1, 2)-submatrix of `A`, and a
    ! third from the (2, 1)-submatrix!
    ! This sounds a little complicated, but it means that we've saved some
    ! memory usage by having two matrices share an object -- the graph g --
    ! rather than duplicate it.
    if (verbose) then
        print *, "o Done creating couplings between (1, 1)- and (2, 2)-"
        print *, "  blocks of A via another random graph g."
        print *, "    Number of references to g:", g%reference_count
    endif



    ! Destroy any heap-allocated objects
    call A%destroy()
    call g%destroy()
    deallocate(g)


! Auxiliary subroutines
contains


!--------------------------------------------------------------------------!
subroutine erdos_renyi_graph(g, m, n, p, symmetric)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(inout) :: g
    integer, intent(in) :: m, n
    real(dp), intent(in) :: p
    logical, intent(in), optional :: symmetric
    ! local variables
    integer :: i, j
    real(dp) :: q
    logical :: sym

    sym = .false.
    if (present(symmetric)) sym = symmetric

    call g%init(m, n)

    do i = 1, m
        if (sym) call g%add_edge(i, i)

        do j = i + 1, n
            call random_number(q)
            if (q < p) then
                call g%add_edge(i, j)
                if (sym) call g%add_edge(j, i)
            endif
        enddo
    enddo

end subroutine erdos_renyi_graph



!--------------------------------------------------------------------------!
subroutine erdos_renyi_matrix(A, g)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(graph_interface), intent(in) :: g
    ! local variables
    integer :: i, j, k, d
    integer, allocatable :: nodes(:)
    real(dp) :: q

    call A%init(g%m, g%n, g)

    d = g%max_degree()
    allocate(nodes(d))

    do i = 1, g%m
        d = g%degree(i)
        call g%get_neighbors(nodes, i)

        do k = 1, d
            j = nodes(k)

            if (j > i) then
                call random_number(q)
                call A%add_value(i, j, -q)
                call A%add_value(i, i, +q)
                call A%add_value(j, i, -q)
                call A%add_value(j, j, +q)
            endif
        enddo
    enddo

end subroutine erdos_renyi_matrix



end program matrix_test_composite

