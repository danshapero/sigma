!==========================================================================!
!==========================================================================!
module ldu_solvers                                                         !
!==========================================================================!
!==========================================================================!
!====     This module contains the definition of LDU solvers. This     ====!
!==== object can be used in more than one way:                         ====!
!====     o A full matrix factorization in which `A = L*D*U`. The      ====!
!====       resulting solver object can then directly solve a linear   ====!
!====       system, to within machine precision and the condition      ====!
!====       number of the underlying matrix.                           ====!
!====     o An incomplete factorization `A ~= L*D*U`. The solver       ====!
!====       is used as a preconditioner for an iterative solver.       ====!
!==== Either way, the representation of the underlying object is the   ====!
!==== same, so we have kept both in the same class.                    ====!
!====     We have implemented a factorization as `A ~= L*D*U`, where   ====!
!==== `L` and `U` are respectively unit lower- and unit upper-         ====!
!==== triangular, and `D` is diagonal. This factorization is most      ====!
!==== convenient, because `U = transpose(L)` when `A` is symmetric.    ====!
!==========================================================================!
!==========================================================================!

use types, only: dp
use linear_operator_interface
use sparse_matrices
use graphs

implicit none




!--------------------------------------------------------------------------!
type, extends(linear_solver) :: sparse_ldu_solver                          !
!--------------------------------------------------------------------------!
    ! Unit triangular factors
    type(csr_matrix) :: L, U

    ! Diagonal entries
    real(dp), allocatable :: D(:)

    ! Permutation of matrix entries
    integer, allocatable :: p(:)

    ! Parameters for solver behavior
    logical :: incomplete = .false.
    logical :: cholesky = .false.
    integer :: level

    logical :: params_set = .false.
contains
    ! Methods required by the linear solver interface
    procedure :: setup => sparse_ldu_setup
    procedure :: linear_solve => ldu_solve
    procedure :: destroy => ldu_destroy

    ! Methods specific to LDU solvers
    procedure :: set_params => ldu_set_params
end type sparse_ldu_solver



contains




!==========================================================================!
!==== Factory method for LDU solvers                                   ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ldu(incomplete,level)                                             !
!--------------------------------------------------------------------------!
    logical, intent(in), optional :: incomplete
    integer, intent(in), optional :: level
    class(linear_solver), pointer :: ldu

    allocate(sparse_ldu_solver::ldu)
    select type(ldu)
        type is(sparse_ldu_solver)
            call ldu%set_params(incomplete, level)
    end select

end function ldu




!==========================================================================!
!==== Core methods for LDU solvers                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_ldu_setup(solver, A)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_ldu_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    ! local variables

    solver%nn = A%nrow

    if (A%ncol /= A%nrow) then
        print *, 'Cannot make an LDU solver for a non-square matrix.'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! Check to make sure `A` is a sparse matrix, not a general linear
    ! operator
    select type(A)
        class is(sparse_matrix_interface)
            if (.not. solver%initialized) then
                ! Build the sparsity pattern for the factorization if we
                ! haven't already done so
                call incomplete_ldu_sparsity_pattern(A, &
                                            & solver%L, solver%U, 0)
                allocate(solver%D(solver%nn))

                solver%initialized = .true.
            endif

            ! Fill the matrix factors
            call sparse_static_pattern_ldu_factorization(A, &
                                            & solver%L, solver%D, solver%U)
    end select


end subroutine sparse_ldu_setup



!--------------------------------------------------------------------------!
subroutine ldu_set_params(solver, incomplete, level)                       !
!--------------------------------------------------------------------------!
    class(sparse_ldu_solver), intent(inout) :: solver
    logical, intent(in), optional :: incomplete
    integer, intent(in), optional :: level

    ! Check if we're performing a full or incomplete factorization;
    ! default to a complete factorization
    !TODO I don't have a direct solver yet, so no matter what make it an
    ! incomplete factorization
    solver%incomplete = .true.

    ! If we're performing an incomplete factorization, find the level of
    ! fill-in allowed
    !TODO I haven't even implemented higher levels of fill, so no matter
    ! what make the level 0
    solver%level = 0

    solver%params_set = .true.

end subroutine ldu_set_params



!--------------------------------------------------------------------------!
subroutine ldu_solve(solver, A, x, b)                                      !
!--------------------------------------------------------------------------!
    class(sparse_ldu_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in)    :: b(:)

    associate(L => solver%L, U => solver%U, D => solver%D)

    x = b
    call lower_triangular_solve(L, x)
    x = x / D
    call upper_triangular_solve(U, x)

    end associate

end subroutine ldu_solve



!--------------------------------------------------------------------------!
subroutine ldu_destroy(solver)                                             !
!--------------------------------------------------------------------------!
    class(sparse_ldu_solver), intent(inout) :: solver

    call solver%L%destroy()
    call solver%U%destroy()

    deallocate(solver%D)

    solver%initialized = .false.

end subroutine ldu_destroy



!==========================================================================!
!==========================================================================!
!==== Auxiliary procedures                                             ====!
!==========================================================================!
!==========================================================================!


!===================================================
!== Forward- and back-solves for triangular systems
!===================================================

!--------------------------------------------------------------------------!
subroutine lower_triangular_solve(M, x)                                    !
!--------------------------------------------------------------------------!
!     Given a strictly lower-triangular matrix `M`, this subroutine solves !
! the linear system `(I+M)x = b`. It is assumed that the right-hand side   !
! vector `b` is stored in the vector `x` to start with.                    !
!     No error-checking is performed in order to make absolutely sure that !
! `M` is strictly lower triangular. If this routine is called on a matrix  !
! which is not strictly lower triangular, it will produce garbage results. !
! However, it is a private routine of this module, to be used only for the !
! matrices constructed herein. If you fall afoul by using this routine     !
! then you must have been very determined to break something.              !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(csr_matrix), intent(in) :: M
    real(dp), intent(inout) :: x(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    do i = 1, M%nrow
        z = x(i)

        do k = M%g%ptr(i), M%g%ptr(i + 1) - 1
            j = M%g%node(k)
            z = z - M%val(k) * x(j)
        enddo

        x(i) = z
    enddo

end subroutine lower_triangular_solve



!--------------------------------------------------------------------------!
subroutine upper_triangular_solve(M, x)                                    !
!--------------------------------------------------------------------------!
!     See comment for previous subroutine.                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(csr_matrix), intent(in) :: M
    real(dp), intent(inout) :: x(:)
    ! local variables
    integer :: i, j, k
    real(dp) :: z

    do i = M%nrow, 1, -1
        z = x(i)

        do k = M%g%ptr(i), M%g%ptr(i + 1) - 1
            j = M%g%node(k)
            z = z - M%val(k) * x(j)
        enddo

        x(i) = z
    enddo

end subroutine upper_triangular_solve




!============================
!== Building the LDU factors
!============================

!--------------------------------------------------------------------------!
subroutine sparse_static_pattern_ldu_factorization(A, L, D, U)             !
!--------------------------------------------------------------------------!
!     This routine takes in a sparse matrix and computes a factorization   !
! into L, D, U, where L is lower triangular, U is upper triangular, and D  !
! is diagonal.                                                             !
!     It is assumed that the sparsity patterns of L and U, namely the      !
! graphs that they point to, have already been created and filled. The     !
! sparsity patterns could correspond to either a complete or incomplete    !
! factorization of the matrix: this subroutine is used the same way in     !
! either case!                                                             !
!     Other subroutines in this module can be used to compute the graphs   !
! for complete or incomplete LDU factorizations.                           !
!     Most algorithms for LDU compute L and U as unit lower triangular     !
! matrices. This algorithm actually computes L-I and U-I, as it is a waste !
! of memory to store n doubles which we know to all be 1.0.                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    type(csr_matrix), intent(inout) :: L, U
    real(dp), intent(inout) :: D(:)
    ! local variables
    integer :: i, j, k, n, nn, dl, du, ind1, ind2
    integer, allocatable :: lneighbors(:), uneighbors(:)
    real(dp) :: Lik, Uki, Uik, Ukj
    ! graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_returned, edges(2, batch_size)
    real(dp) :: entries(batch_size)

    nn = A%nrow

    call L%zero()
    call U%zero()
    D = 0.0_dp

    ! Copy A into L, D, U
    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_entries(edges, entries, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            if (i > j) then
                call L%set_value(i, j, entries(k))
            elseif (j > i) then
                call U%set_value(i, j, entries(k))
            else
                D(i) = entries(k)
            endif
        enddo
    enddo

    dl = L%g%max_degree()
    du = U%g%max_degree()

    allocate(lneighbors(dl), uneighbors(du))

    do i = 1, nn
        ! Get the degree and neighbors of node `i` in `L` and `U`
        call L%g%get_neighbors(lneighbors, i)
        dl = L%g%degree(i)

        call U%g%get_neighbors(uneighbors, i)
        du = U%g%degree(i)

        do ind1 = 1, dl
            k = lneighbors(ind1)

            Lik = L%get_value(i, k)
            Uki = U%get_value(k, i)

            call L%set_value(i, k, Lik / D(k))
            Lik = Lik / D(k)

            ! Loop through all neighbors `j` of vertex `i` in `L` with
            ! `j > k`...
            do ind2 = 1, dl
                j = lneighbors(ind2)

                if (j > k) then
                    !... and add up their contribution to `L`.
                    Ukj = U%get_value(k, j)
                    call L%add_value(i, j, -Lik * D(k) * Ukj)
                endif
            enddo

            ! Modify the diagonal entry
            D(i) = D(i) - Lik * D(k) * Uki

            ! Now loop through all neighbors `j` of vertex `i` in `U`
            do ind2 = 1, du
                j = uneighbors(ind2)
                Ukj = U%get_value(k, j)
                call U%add_value(i, j, -Lik * D(k) * Ukj)
            enddo
        enddo   ! End of loop over neighbors `k` of `i`

        ! Finally, loop through row `i` in `U` and scale all the entries by
        ! the diagonal
        do ind2 = 1, du
            k = uneighbors(ind2)
            Uik = U%get_value(i, k)
            call U%set_value(i, k, Uik / D(i))
        enddo

    enddo   ! End of loop over all rows `i`


    deallocate(lneighbors, uneighbors)

end subroutine sparse_static_pattern_ldu_factorization




!==================================
!== Computing the sparsity pattern
!==================================

!--------------------------------------------------------------------------!
subroutine incomplete_ldu_sparsity_pattern(A, L, U, level)                 !
!--------------------------------------------------------------------------!
!     This subroutine builds the sparse matrix factors `L`, `U` for an     !
! incomplete LDU factorization of the matrix `A`. The argument `level`     !
! specifies the level of fill allowed in the factorization.                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    type(csr_matrix), intent(inout) :: L, U
    integer, intent(in) :: level
    ! local variables
    integer :: i, j, k, nn
    type(ll_graph) :: gl, gu
    ! graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: edges(2, batch_size), num_returned

    if (level > 0) then
        print *, 'ILDU(k) for k > 0 not implemented yet! Sorry pal.'
        call exit(1)
    endif

    nn = A%nrow

    call gl%init(nn)
    call gu%init(nn)

    cursor = A%make_cursor()
    do while (.not. cursor%done())
        call A%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            if (i > j) call gl%add_edge(i, j)
            if (j > i) call gu%add_edge(i, j)
        enddo
    enddo

    call L%init(nn, nn, gl)
    call U%init(nn, nn, gu)

end subroutine incomplete_ldu_sparsity_pattern



end module ldu_solvers

