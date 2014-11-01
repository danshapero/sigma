!==========================================================================!
!==========================================================================!
module sparse_matrix_algebra                                               !
!==========================================================================!
!==========================================================================!

use sparse_matrix_interfaces
use cs_matrices
use ll_graphs

implicit none

private
public :: sparse_matrix_sum, sparse_matrix_product, PtAP, RARt

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
! This is a driver routine for computing the product                       !
!     A = B * C                                                            !
! of two sparse matrices; it selects the proper implementation to use      !
! based on whether the factors are stored in row / column major format.    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    type(csr_matrix) :: CC
    type(ll_graph) :: g

    if (B%ncol /= C%nrow) then
        print *, 'Dimensions of matrices for sparse matrix product'
        print *, 'A = B * C are not consistent.'
        print *, 'Column dimension of B:', B%ncol
        print *, 'Row dimension of C:   ', C%nrow
        call exit(1)
    endif

    call A%set_dimensions(B%nrow, C%ncol)

    if (C%is_get_row_fast()) then
        call sparse_matrix_product_C(A, B, C)
    elseif (B%is_get_column_fast()) then
        call sparse_matrix_product_B(A, B, C)
    else
        call CC%set_dimensions(B%nrow, C%ncol)
        call CC%copy_matrix(C)
        call sparse_matrix_product_C(A, B, CC)
    endif

end subroutine sparse_matrix_product



!--------------------------------------------------------------------------!
subroutine sparse_matrix_product_B(A, B, C)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    type(ll_graph) :: g

    call matrix_product_graph(g, B, C)
    call A%copy_graph(g)
    call matrix_product_fill_entries(A, B, C)

contains

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine matrix_product_graph(g, B, C)                               !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
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
        d = B%get_max_column_degree()
        allocate(nodes(d), slice(d))

        ! Iterate through all the edges of `C`
        cursor = C%make_cursor()
        do while(.not. cursor%done())
            call C%get_edges(edges, cursor, batch_size, num_returned)

            do l = 1, num_returned
                ! For each edge (k, j) in `C`
                k = edges(1, l)
                j = edges(2, l)

                ! Find all edges (i, k) in `B` and add (i, j) to `g`
                call B%get_column(nodes, slice, k)
                d = B%get_column_degree(k)

                do m = 1, d
                    i = nodes(m)
                    call g%add_edge(i, j)
                enddo
            enddo
        enddo

        deallocate(nodes, slice)

    end subroutine matrix_product_graph


    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine matrix_product_fill_entries(A, B, C)                        !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
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

        ! Allocate arrays for getting slices from `B`
        d = B%get_max_column_degree()
        allocate(nodes(d), slice(d))

        ! Iterate through all the edges of `C`
        cursor = C%make_cursor()
        do while(.not. cursor%done())
            call C%get_entries(edges, vals, cursor, batch_size, num_returned)

            do l = 1, num_returned
                ! For each edge (k, j) in `C`
                k = edges(1, l)
                j = edges(2, l)
                Ckj = vals(l)

                ! Find all edges (i, k) in `B` and add B(i, k) * C(k, j)
                ! to A(i, j)
                call B%get_column(nodes, slice, k)
                d = B%get_column_degree(k)

                do m = 1, d
                    i = nodes(m)
                    Bik = slice(m)

                    call A%add_value(i, j, Bik * Ckj)
                enddo
            enddo
        enddo

    end subroutine matrix_product_fill_entries

end subroutine sparse_matrix_product_B



!--------------------------------------------------------------------------!
subroutine sparse_matrix_product_C(A, B, C)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: A
    class(sparse_matrix_interface), intent(in) :: B, C
    ! local variables
    type(ll_graph) :: g

    call matrix_product_graph(g, B, C)
    call A%copy_graph(g)
    call matrix_product_fill_entries(A, B, C)

contains

    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine matrix_product_graph(g, B, C)                               !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
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
        d = C%get_max_row_degree()
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

    end subroutine matrix_product_graph


    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    subroutine matrix_product_fill_entries(A, B, C)                        !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
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
        d = C%get_max_row_degree()
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

    end subroutine

end subroutine sparse_matrix_product_C



!--------------------------------------------------------------------------!
subroutine PtAP(B, A, P)                                                   !
!--------------------------------------------------------------------------!
! This procedure computes the congruent product                            !
!     B = P^t * A * P                                                      !
! of the sparse matrices A, P. This operation is common in multigrid and   !
! domain-decomposition methods. It is assumed that the matrix `P` is in a  !
! format for which getting a matrix row can done quickly, otherwise this   !
! operation will be very slow.                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(inout) :: B
    class(sparse_matrix_interface), intent(in) :: A, P
    ! local variables
    integer :: i, j, k, l, m, d
    real(dp) :: Akl, Pki, Plj
    ! variables for getting row slices from P
    integer :: n1, n2, d1, d2
    integer, allocatable :: nodes1(:), nodes2(:)
    real(dp), allocatable :: slice1(:), slice2(:)
    ! variables for iterating through the values of A
    type(graph_edge_cursor) :: cursor
    integer :: num_returned, edges(2, batch_size)
    real(dp) :: vals(batch_size)
    ! graph for the connectivity structure of B
    type(ll_graph) :: g

    ! Check that all matrices have the right dimensions
    if (A%nrow /= A%ncol) then
        print *, "Cannot perform operation B = P^t * A * P for a non-"
        print *, "square matrix A!"
        call exit(1)
    endif

    if (A%ncol /= P%nrow) then
        print *, "Dimension mismatch in operation B = P^t * A * P."
        print *, "Dimension of A:", A%nrow, A%ncol
        print *, "Dimension of P:", P%nrow, P%ncol
        call exit(1)
    endif

    d = P%get_max_row_degree()
    allocate(nodes1(d), nodes2(d), slice1(d), slice2(d))

    call g%init(P%ncol, P%ncol)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! First, build the connectivity graph of B by iterating over all the   !
    ! edges of `A` and using the relation                                  !
    !     B_ij = sum_{k, l} P_ki * A_kl * P_lj                             !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_edges(edges, cursor, batch_size, num_returned)

        do m = 1, num_returned
            k = edges(1, m)
            l = edges(2, m)

            d1 = P%get_row_degree(k)
            call P%get_row(nodes1, slice1, k)

            d2 = P%get_row_degree(l)
            call P%get_row(nodes2, slice2, l)

            do n1 = 1, d1
                i = nodes1(n1)

                do n2 = 1, d2
                    j = nodes2(n2)
                    call g%add_edge(i, j)
                enddo
            enddo
        enddo
    enddo

    call B%init(P%ncol, P%ncol, g)
    call B%zero()

    call g%destroy()


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! Fill the entries of B                                                !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_entries(edges, vals, cursor, batch_size, num_returned)

        do m = 1, num_returned
            k = edges(1, m)
            l = edges(2, m)
            Akl = vals(m)

            d1 = P%get_row_degree(k)
            call P%get_row(nodes1, slice1, k)

            d2 = P%get_row_degree(l)
            call P%get_row(nodes2, slice2, l)

            do n1 = 1, d1
                i = nodes1(n1)
                Pki = slice1(n1)

                do n2 = 1, d2
                    j = nodes2(n2)
                    Plj = slice2(n2)

                    call B%add(i, j, Pki * Akl * Plj)
                enddo
            enddo
        enddo
    enddo

end subroutine PtAP



!--------------------------------------------------------------------------!
subroutine RARt(B, A, R)                                                   !
!--------------------------------------------------------------------------!
! This procedure computes the congruent product                            !
!     B = R * A * R^t                                                      !
! of the sparse matrices A, R. This operation is common in multigrid and   !
! domain-decomposition methods. It is assumed that the matrix `R` is in a  !
! format for which getting a matrix column can done quickly, otherwise     !
! this operation will be very slow.                                        !
!--------------------------------------------------------------------------!
    class(sparse_matrix_interface), intent(inout) :: B
    class(sparse_matrix_interface), intent(in) :: A, R
    ! local variables
    integer :: i, j, k, l, m, d
    real(dp) :: Akl, Rik, Rjl
    ! variables for getting row slices from P
    integer :: n1, n2, d1, d2
    integer, allocatable :: nodes1(:), nodes2(:)
    real(dp), allocatable :: slice1(:), slice2(:)
    ! variables for iterating through the values of A
    type(graph_edge_cursor) :: cursor
    integer :: num_returned, edges(2, batch_size)
    real(dp) :: vals(batch_size)
    ! graph for the connectivity structure of B
    type(ll_graph) :: g

    ! Check that all matrices have the right dimensions
    if (A%nrow /= A%ncol) then
        print *, "Cannot perform operation B = P^t * A * P for a non-"
        print *, "square matrix A!"
        call exit(1)
    endif

    if (R%ncol /= A%nrow) then
        print *, "Dimension mismatch in operation B = R * A * R^t."
        print *, "Dimension of A:", A%nrow, A%ncol
        print *, "Dimension of P:", R%nrow, R%ncol
        call exit(1)
    endif

    d = R%get_max_column_degree()
    allocate(nodes1(d), nodes2(d), slice1(d), slice2(d))

    call g%init(R%nrow, R%nrow)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! First, build the connectivity graph of B by iterating over all the   !
    ! edges of `A` and using the relation                                  !
    !     B_ij = sum_{k, l} R_ik * A_kl * R_jl                             !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_edges(edges, cursor, batch_size, num_returned)

        do m = 1, num_returned
            k = edges(1, m)
            l = edges(2, m)

            d1 = R%get_column_degree(k)
            call R%get_column(nodes1, slice1, k)

            d2 = R%get_column_degree(l)
            call R%get_column(nodes2, slice2, l)

            do n1 = 1, d1
                i = nodes1(n1)

                do n2 = 1, d2
                    j = nodes2(n2)
                    call g%add_edge(i, j)
                enddo
            enddo
        enddo
    enddo

    call B%init(R%nrow, R%nrow, g)
    call B%zero()

    call g%destroy()


    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
    ! Fill the entries of B                                                !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_entries(edges, vals, cursor, batch_size, num_returned)

        do m = 1, num_returned
            k = edges(1, m)
            l = edges(2, m)
            Akl = vals(m)

            d1 = R%get_column_degree(k)
            call R%get_column(nodes1, slice1, k)

            d2 = R%get_column_degree(l)
            call R%get_column(nodes2, slice2, l)

            do n1 = 1, d1
                i = nodes1(n1)
                Rik = slice1(n1)

                do n2 = 1, d2
                    j = nodes2(n2)
                    Rjl = slice2(n2)

                    call B%add(i, j, Rik * Akl * Rjl)
                enddo
            enddo
        enddo
    enddo

end subroutine RARt



end module sparse_matrix_algebra

