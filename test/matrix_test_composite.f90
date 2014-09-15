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
    type(ll_graph) :: g

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
    nn = 768
    p = log(1.0_dp * nn) / log(2.0_dp) / nn

    ! Make an Erdos-Renyi graph
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        do j = i + 1, nn
            call random_number(q)

            if (q < p) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    d = g%max_degree()

    if (verbose) then
        print *, "o Done generating Erdos-Renyi graph."
        print *, "    Number of vertices: ", nn
        print *, "    Number of edges:    ", g%ne
        print *, "    Maximum degree:     ", d
    endif


    allocate(nodes(d))

    ! Make a Laplacian matrix on the graph with random weights
    call choose_matrix_type(C, "csr")
    call C%init(nn, nn, g)
    call C%zero()

    do i = 1, nn
        d = g%degree(i)
        call g%get_neighbors(nodes, i)

        do k = 1, d
            j = nodes(k)
            if (j > i) then
                call random_number(q)
                call C%add_value(i, i, +q)
                call C%add_value(j, j, +q)
                call C%add_value(i, j, -q)
                call C%add_value(j, i, -q)
            endif
        enddo
    enddo

    if (verbose) then
        print *, "o Done generating random weighted Laplacian matrix."
    endif


    ! Set the (1, 1)-submatrix of `A` to be `C`
    call A%set_submatrix(1, 1, C)



    call A%destroy()
    call C%destroy()
    deallocate(C)

end program matrix_test_composite

