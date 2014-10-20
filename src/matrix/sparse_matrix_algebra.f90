!==========================================================================!
!==========================================================================!
module sparse_matrix_algebra                                               !
!==========================================================================!
!==========================================================================!

use sparse_matrix_interfaces
use ll_graphs

implicit none

private
public :: sparse_matrix_sum, sparse_matrix_product

contains



!==========================================================================!
!==== Routines for matrix sums                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_matrix_sum(A, B, C)                                      !
!--------------------------------------------------------------------------!
!     This routine computes the sum `A` of the sparse matrices `B`, `C`.   !
! The actual work is passed off to two procedures (see below) that compute !
! the connectivity structure of the matrix sum, then fill the entries.     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in)    :: B, C
    ! local variables
    type(ll_graph) :: g

    if (B%nrow /= C%nrow .or. B%ncol /= C%ncol) then
        print *, 'Dimensions of matrices for sparse matrix sum'
        print *, 'A = B+C are not consistent.'
        print *, 'Dimensions of B:', B%nrow, B%ncol
        print *, 'Dimensions of C:', C%nrow, C%ncol
        print *, 'Terminating.'
        call exit(1)
    endif

    call A%set_dimensions(B%nrow, B%ncol)
    call sparse_matrix_sum_graph(g, B, C)
    call A%copy_graph(g)
    call sparse_matrix_sum_fill_entries(A, B, C)

end subroutine sparse_matrix_sum



!--------------------------------------------------------------------------!
subroutine sparse_matrix_sum_graph(g, B, C)                                !
!--------------------------------------------------------------------------!
!     This routine computes the graph `g` describing the underlying        !
! connectivity structure of the sparse matrix sum `B + C`. It is a helper  !
! routine for other procedures which compute the entries of the sum.       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(inout) :: g
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    integer :: i, j, k
    integer :: num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    call g%init(B%nrow, B%ncol)


    ! Add all the edges from `B` into `g`
    cursor = B%make_cursor()
    do while(.not. cursor%done())
        call B%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            call g%add_edge(i, j)
        enddo
    enddo


    ! Now add all the edges from `C` into `g`
    cursor = C%make_cursor()
    do while(.not. cursor%done())
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
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in)    :: B, C
    ! local variables
    integer :: i, j, k
    real(dp) :: z
    integer :: num_returned, edges(2, batch_size)
    real(dp) :: entries(batch_size)
    type(graph_edge_cursor) :: cursor

    ! Add up the contributions to `A` from `B`
    cursor = B%make_cursor()
    do while(.not. cursor%done())
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
    do while (.not. cursor%done())
        call C%get_entries(edges, entries, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)
            z = entries(k)

            call A%add_value(i, j, z)
        enddo
    enddo

end subroutine sparse_matrix_sum_fill_entries



!==========================================================================!
!==== Routines for matrix products                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_matrix_product(A, B, C)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    type(ll_graph) :: g

    if (B%ncol /= C%nrow) then
        print *, 'Dimensions of matrices for sparse matrix product'
        print *, 'A = B * C are not consistent.'
        print *, 'Column dimension of B:', B%ncol
        print *, 'Row dimension of C:   ', C%nrow
        call exit(1)
    endif

    call A%set_dimensions(B%nrow, C%ncol)
    call sparse_matrix_product_graph(g, B, C)
    call A%copy_graph(g)
    call sparse_matrix_product_fill_entries(A, B, C)

end subroutine sparse_matrix_product



!--------------------------------------------------------------------------!
subroutine sparse_matrix_product_graph(g, B, C)                            !
!--------------------------------------------------------------------------!
    class(graph_interface), intent(inout) :: g
    class(sparse_matrix_interface), intent(in) :: B, C

    if (C%is_get_row_fast()) then
        call product_graph_C(g, B, C)
    elseif (B%is_get_column_fast()) then
        call product_graph_B(g, B, C)
    else
        ! In this latter case, the first factor does not have fast column
        ! access and the second factor does not have fast row access --
        ! this is *very* bad and we should probably copy one of the matrices
        ! to either CSR or CSC and compute the product using the copy.
        call product_graph_C(g, B, C)
    endif

contains
    ! These nested procedures do the actual work to compute the graph
    ! representing the connectivity structure of the product matrix.

    !-----------------------------------!
    subroutine product_graph_C(g, B, C) !
        ! input/output variables
        class(graph_interface), intent(inout) :: g
        class(sparse_matrix_interface), intent(in) :: B, C
        ! local variables
        integer :: i, j, k, l, m, d
        ! variables for getting matrix slices
        integer, allocatable :: nodes(:)
        real(dp), allocatable :: slice(:)
        ! variables for iterating through matrix entries
        integer :: num_returned, edges(2, batch_size)
        type(graph_edge_cursor) :: cursor

        call g%init(B%nrow, C%ncol)

        ! Allocate arrays for getting slices from `C`
    !    d = C%max_row_degree() ! Haha this doesn't even exist dumbass
        d = 0
        do i = 1, C%nrow
            d = max(d, C%get_row_degree(i))
        enddo
        allocate(nodes(d), slice(d))


        ! Iterate through all the edges of `B`
        cursor = B%make_cursor()
        do while(.not. cursor%done())
            call B%get_edges(edges, cursor, batch_size, num_returned)

            do l = 1, num_returned
                ! For each edge (i, k) in `B`
                i = edges(1, l)
                k = edges(2, l)

                ! Find all edges (k, j) in `C` and add (i, j) to `g`
                call C%get_row(nodes, slice, k)
                d = C%get_row_degree(k)

                do m = 1, d
                    j = nodes(m)
                    call g%add_edge(i, j)
                enddo
            enddo
        enddo

        deallocate(nodes, slice)

    end subroutine product_graph_C

    !-----------------------------------!
    subroutine product_graph_B(g, B, C) !
        ! input/output variables
        class(graph_interface), intent(inout) :: g
        class(sparse_matrix_interface), intent(in) :: B, C
        ! local variables
        integer :: i, j, k, l, m, d
        ! variables for getting matrix slices
        integer, allocatable :: nodes(:)
        real(dp), allocatable :: slice(:)
        ! variables for iterating through matrix entries
        integer :: num_returned, edges(2, batch_size)
        type(graph_edge_cursor) :: cursor

        call g%init(B%nrow, C%ncol)

        ! Allocate arrays for getting slices from `B`
    !    d = C%max_column_degree() ! Haha this doesn't even exist dumbass
        d = 0
        do j = 1, B%ncol
            d = max(d, B%get_column_degree(j))
        enddo
        allocate(nodes(d), slice(d))

        ! Iterate through all the edges of `B`
        cursor = C%make_cursor()
        do while(.not. cursor%done())
            call C%get_edges(edges, cursor, batch_size, num_returned)

            do l = 1, num_returned
                ! For each edge (i, k) in `B`
                k = edges(1, l)
                j = edges(2, l)

                ! Find all edges (k, j) in `C` and add (i, j) to `g`
                call B%get_column(nodes, slice, k)
                d = B%get_column_degree(k)

                do m = 1, d
                    i = nodes(m)
                    call g%add_edge(i, j)
                enddo
            enddo
        enddo

        deallocate(nodes, slice)

    end subroutine product_graph_B

end subroutine sparse_matrix_product_graph



!--------------------------------------------------------------------------!
subroutine sparse_matrix_product_fill_entries(A, B, C)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    integer :: i, j, k, l, m, d
    real(dp) :: Bik, Ckj
    ! variables for getting matrix slices
    integer, allocatable :: nodes(:)
    real(dp), allocatable :: slice(:)
    ! variables for iterating through matrix entries
    integer :: num_returned, edges(2, batch_size)
    real(dp) :: vals(batch_size)
    type(graph_edge_cursor) :: cursor

    call A%zero()

    ! Allocate arrays for getting slices from `C`
    !d = C%max_row_degree()
    d = 0
    do i = 1, C%nrow
        d = max(d, C%get_row_degree(i))
    enddo
    allocate(nodes(d), slice(d))


    ! Iterate through all the edges of `B`
    cursor = B%make_cursor()
    do while(.not. cursor%done())
        call B%get_entries(edges, vals, cursor, batch_size, num_returned)

        do l = 1, num_returned
            ! For each edge (i, k) in `B`
            i = edges(1, l)
            k = edges(2, l)
            Bik = vals(l)

            ! Find all edges (k, j) in `C` and add B(i, k) * C(k, j) 
            ! to A(i, j)
            call C%get_row(nodes, slice, k)
            d = C%get_row_degree(k)
            do m = 1, d
                j = nodes(m)
                Ckj = slice(m)

                call A%add_value(i, j, Bik * Ckj)
            enddo
        enddo
    enddo    

end subroutine sparse_matrix_product_fill_entries



end module sparse_matrix_algebra

