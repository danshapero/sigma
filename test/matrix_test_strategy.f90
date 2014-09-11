!--------------------------------------------------------------------------!
program matrix_test_strategy                                               !
!--------------------------------------------------------------------------!
!     This program tests using the type sparse_matrix as the container of  !
! a strategy. The abstract strategy is                                     !
!     sparse_matrix_interface                                              !
! and the concrete strategy is one of the sparse matrix storage formats,   !
! such as CSR or ellpack.                                                  !
!     Using this object allows the user to transparently change the        !
! storage format of a sparse matrix without having to explicitly change    !
! the dynamic type of the object. Instead, that occurs behind the scenes.  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    type(ll_graph) :: g

    ! sparse matrix objects
    type(sparse_matrix) :: A

    ! vectors
    real(dp), allocatable :: x(:), y(:)

    ! integer indices
    integer :: i, j, k, d, nn

    ! permutation
    integer, allocatable :: p(:)

    ! variables for getting matrix rows / columns
    integer :: row_degree, col_degree
    integer, allocatable :: nodes(:)
    real(dp), allocatable :: slice(:)

    ! random numbers
    real(dp) :: c, w, z

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
    c = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make a random reference sparse matrix                                !
    !----------------------------------------------------------------------!
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        do j = i + 1, nn
            call random_number(z)

            if (z < c) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    call A%set_matrix_type("csr")
    call A%set_dimensions(nn, nn)
    call A%copy_graph_structure(g)


    call g%destroy()
    call A%destroy()


end program matrix_test_strategy

