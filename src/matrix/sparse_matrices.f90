!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!

use graphs
use sparse_matrix_interface
use default_sparse_matrix_kernels
use default_matrices
use cs_matrices
use ellpack_matrices

implicit none

interface sparse_matrix
    module procedure sparse_matrix_factory
end interface


contains



!--------------------------------------------------------------------------!
function sparse_matrix_factory(nrow, ncol, g, orientation) result(A)       !
!--------------------------------------------------------------------------!
    integer, intent(in) :: nrow, ncol
    class(graph), pointer, intent(in) :: g
    character(len=3), intent(in) :: orientation
    class(sparse_matrix), pointer :: A

    select type(g)
        class is(cs_graph)
            A => cs_matrix(nrow, ncol, g, orientation)
        class is(ellpack_graph)
            A => ellpack_matrix(nrow, ncol, g, orientation)
        class default
            A => default_matrix(nrow, ncol, g, orientation)
    end select

end function sparse_matrix_factory




!==========================================================================!
!==== Algebraic operations on sparse matrices                          ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_matrix_sum_graph(g, B, C)                                !
!--------------------------------------------------------------------------!
!     This routines computes the graph `g` describing the underlying       !
! connectivity structure of the sparse matrix sum `B + C`. It is a helper  !
! routine for other procedures which compute the entries of the sum.       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    class(sparse_matrix), intent(in) :: B, C
    ! local variables
    integer :: i, j, k
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    call g%init(B%nrow, B%ncol)


    ! Add all the edges from `B` into `g`
    cursor = B%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call B%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            call g%add_edge(i, j)
        enddo
    enddo


    ! Now add all the edges from `C` into `g`
    cursor = C%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call C%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            call g%add_edge(i, j)
        enddo
    enddo

end subroutine sparse_matrix_sum_graph



!--------------------------------------------------------------------------!
subroutine sparse_matrix_sum_fill_entries(A, B, C)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in)    :: B, C
    ! local variables
    integer :: i, j, k
    real(dp) :: z
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    real(dp) :: entries(batch_size)
    type(graph_edge_cursor) :: cursor

    ! Add up the contributions to `A` from `B`
    cursor = B%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call B%get_entries(edges, entries, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)
            z = entries(k)

            call A%add_value(i, j, z)
        enddo
    enddo

    ! Add up the contributions to `A` from `C`
    cursor = C%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call C%get_entries(edges, entries, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)
            z = entries(k)

            call A%add_value(i, j, z)
        enddo
    enddo

end subroutine sparse_matrix_sum_fill_entries





end module sparse_matrices
