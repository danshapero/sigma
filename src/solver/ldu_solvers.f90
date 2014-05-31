module ldu_solvers

use types, only: dp
use linear_operator_interface
use sparse_matrices
use graph_interface
use cs_graphs

implicit none




!--------------------------------------------------------------------------!
type, extends(linear_solver) :: sparse_ldu_solver                          !
!--------------------------------------------------------------------------!
    ! Triangular factors L and U, stored in one matrix
    !TODO: make these attributes private once testing is complete
    !type(sparse_matrix), private :: L, U
    type(sparse_matrix) :: L, U

    ! Diagonal entries D
    !TODO: make these attributes private once testing is complete
    !real(dp), allocatable, private :: D(:)
    real(dp), allocatable :: D(:)

    ! Scratch vector for intermediate computations
    real(dp), allocatable :: z(:)

    ! Permutation of matrix entries in case of pivoting
    integer, allocatable :: p(:)

    ! Parameters determining behavior of the solver
    logical :: incomplete = .false.
    logical :: cholesky = .false.
    logical :: params_set = .false.
    integer :: level
contains
    ! Methods required by the linear solver interface
    procedure :: setup => sparse_ldu_setup
    procedure :: linear_solve => ldu_solve
    procedure :: destroy => ldu_destroy

    ! Methods specific to LDU solvers
    procedure :: set_params => ldu_set_params
end type sparse_ldu_solver



!TODO: make these methods private once testing is complete
!--------------------------------------------------------------------------!
!private :: lower_triangular_solve, upper_triangular_solve                  !
!private :: sparse_static_pattern_ldu_factorization                         !
!private :: incomplete_ldu_sparsity_pattern                                 !
!--------------------------------------------------------------------------!



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
            call ldu%set_params(incomplete,level)
    end select

end function ldu




!==========================================================================!
!==== Core methods for LDU solvers                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_ldu_setup(solver,A)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_ldu_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    ! local variables
    class(graph), pointer :: g, gL, gU

    solver%nn = A%nrow

    if (A%ncol/=A%nrow) then
        print *, 'Cannot make an LDU solver for a non-square matrix.'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! First, check to make sure A is a sparse matrix, and not block-sparse
    ! or a general linear operator
    select type(A)
        type is(sparse_matrix)
            if (.not.solver%initialized) then
                g => A%g

                ! Build the sparsity pattern for the factorization
                call incomplete_ldu_sparsity_pattern( gL, gU, g, 0)

                ! Initialize the L, U factors and allocate space for the
                ! diagonal elements D
                call solver%L%init(solver%nn,solver%nn,'row',gL)
                call solver%U%init(solver%nn,solver%nn,'row',gU)
                allocate(solver%D(solver%nn))

                solver%initialized = .true.
            endif

            ! Fill the L, D, U factors
            call sparse_static_pattern_ldu_factorization(A, &
                & solver%L, solver%D, solver%U)
    end select

end subroutine sparse_ldu_setup



!--------------------------------------------------------------------------!
subroutine ldu_set_params(solver,incomplete,level)                         !
!--------------------------------------------------------------------------!
    class(sparse_ldu_solver), intent(inout) :: solver
    logical, intent(in), optional :: incomplete
    integer, intent(in), optional :: level

    ! Check if we're performing a full or incomplete factorization;
    ! default to a complete factorization
    if (present(incomplete)) then
        solver%incomplete = incomplete
    else
        solver%incomplete = .false.
    endif

    ! If we're performing an incomplete factorization, find the level of
    ! fill-in allowed
    if (solver%incomplete) then
        if (present(level)) then
            solver%level = level
        else
            solver%level = 0
        endif
    endif

    solver%params_set = .true.

end subroutine



!--------------------------------------------------------------------------!
subroutine ldu_solve(solver,A,x,b)                                         !
!--------------------------------------------------------------------------!
    class(sparse_ldu_solver), intent(inout) :: solver
    class(linear_operator), intent(in) :: A
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in)    :: b(:)

    associate(L => solver%L, U => solver%U, D => solver%D)

    x = b
    call lower_triangular_solve(L,x)
    x = x/D
    call upper_triangular_solve(U,x)

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
subroutine lower_triangular_solve(M,x)                                     !
!--------------------------------------------------------------------------!
!     Given a strictly lower-triangular matrix M, this subroutine solves   !
! the linear system (I+M)x = b. It is assumed that the right-hand side     !
! vector b is stored in the vector x to start with.                        !
!     No error-checking is performed in order to make absolutely sure that !
! M is strictly lower triangular. If this routine is called on a matrix    !
! which is not strictly lower triangular, it will produce garbage results. !
! However, it is a private routine of this module, to be used only for the !
! matrices constructed herein. If you fall afoul by using this routine     !
! then you must have been very determined to break something.              !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(sparse_matrix), intent(in) :: M
    real(dp), intent(inout) :: x(:)
    ! local variables
    integer :: i, j, k, d, neighbors(M%g%max_degree)
    real(dp) :: vals(M%g%max_degree)

    do i=1,M%nrow
        call M%get_row(neighbors,vals,i)

        d = M%g%degree(i)
        do k=1,d
            j = neighbors(k)
            x(i) = x(i)-vals(k)*x(j)
        enddo
    enddo

end subroutine lower_triangular_solve



!--------------------------------------------------------------------------!
subroutine upper_triangular_solve(M,x)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(sparse_matrix), intent(in) :: M
    real(dp), intent(inout) :: x(:)
    ! local variables
    integer :: i, j, k, d, neighbors(M%g%max_degree)
    real(dp) :: vals(M%g%max_degree)

    do i=M%nrow,1,-1
        call M%get_row(neighbors,vals,i)

        d = M%g%degree(i)
        do k=1,d
            j = neighbors(k)
            x(i) = x(i)-vals(k)*x(j)
        enddo
    enddo

end subroutine upper_triangular_solve




!============================
!== Building the LDU factors
!============================

!--------------------------------------------------------------------------!
subroutine sparse_static_pattern_ldu_factorization(A,L,D,U)                !
!--------------------------------------------------------------------------!
!     This routine takes in a sparse matrix A and computes a factorization !
! into L, D, U, where L is lower triangular, U is upper triangular, and D  !
! is diagonal.                                                             !
!     It is assumed that the sparsity patterns of L and U, namely the      !
! graphs that they point to, have already been created and filled. The     !
! sparsity patterns could correspond to either a complete or incomplete    !
! factorization of the matrix: this subroutine is used the same way in     !
! either case!                                                             !
!     Other subroutines in this module can be used to compute the graphs   !
! for complete or incomplete LDU factorizations. For the purposes of this  !
! algorithm, it is assumed that both L and U use an underlying graph data  !
! structure for which finding all the neighbors of a given vertex is fast, !
! and that they are both stored in row-major order, not column-major. This !
! algorithm uses the IKJ variant of Gaussian elimination, which is best    !
! employed when row operations are fast.                                   !
!     Most algorithms for LDU compute L and U as unit lower triangular     !
! matrices. This algorithm actually computes L-I and U-I, as it is a waste !
! of memory to store n doubles which we know to all be 1.0.                !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(sparse_matrix), intent(in) :: A
    type(sparse_matrix), intent(inout) :: L, U
    real(dp), intent(inout) :: D(:)
    ! local variables
    integer :: i,j,k,n,nn,dl,du,ind1,ind2
    integer, allocatable :: lneighbors(:),uneighbors(:)
    real(dp) :: Lik, Uki, Uik, Ukj
    ! graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64), order(2)

    nn = A%nrow

    call L%zero()
    call U%zero()
    D = 0.0_dp

    ! Copy A into L, D and U
    order = A%order
    cursor = A%g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1
    do n=1,num_blocks
        call A%g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            if (i/=0 .and. j/=0) then
                if (i>j) then
                    call L%set_value(i,j,A%val(64*(n-1)+k))
                elseif(j>i) then
                    call U%set_value(i,j,A%val(64*(n-1)+k))
                else
                    D(i) = A%val(64*(n-1)+k)
                endif
            endif
        enddo
    enddo

    ! Allocate a neighbors array
    allocate(lneighbors(L%g%max_degree), uneighbors(U%g%max_degree))

    ! Compute the factorization using the IKJ-variant of Gaussian
    ! elimination
    do i=1,nn
        ! Get the degree and neighbors of node i in L
        call L%g%get_neighbors(lneighbors,i)
        dl = L%g%degree(i)

        ! Get the degree and neighbors of node i in U
        call U%g%get_neighbors(uneighbors,i)
        du = U%g%degree(i)

        ! Loop through all the neighbors k of vertex i in the matrix L
        do ind1=1,dl
            k = lneighbors(ind1)

            ! Get the (i,k) entry of L
            Lik = L%get_value(i,k)
            Uki = U%get_value(k,i)

            ! Divide the L(i,k) by D(k)
            call L%set_value(i,k,Lik/D(k))
            Lik = Lik/D(k)

            ! Now loop through all neighbors j of vertex i in the matrix L
            ! such that j > k
            do ind2=1,dl
                j = lneighbors(ind2)

                ! If j > k,
                if (j>k) then
                    ! Set the value of LU(k,j)
                    Ukj = U%get_value(k,j)
                    call L%add_value(i,j,-Lik*D(k)*Ukj)
                endif
            enddo

            ! Adjust the diagonal
            D(i) = D(i)-Lik*D(k)*Uki

            ! And finally loop through all the neighbors j of vertex i in
            ! the matrix U
            do ind2=1,du
                j = uneighbors(ind2)
                Ukj = U%get_value(k,j)
                call U%add_value(i,j,-Lik*D(k)*Ukj)
            enddo
        enddo

        ! Loop through all the neighbors k of vertex i in the matrix U
        do ind2=1,du
            ! Scale them by the diagonal
            k = uneighbors(ind2)
            Uik = U%get_value(i,k)
            call U%set_value(i,k,Uik/D(i))
        enddo
    enddo

end subroutine sparse_static_pattern_ldu_factorization




!==================================
!== Computing the sparsity pattern
!==================================

!--------------------------------------------------------------------------!
subroutine incomplete_ldu_sparsity_pattern(gL,gU,g,level,trans)            !
!--------------------------------------------------------------------------!
!     This subroutine builds the graphs gL, gU for the factors L, U of an  !
! incomplete LDU factorization of the matrix A whose connectivity is       !
! represented by the graph g.                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), pointer, intent(out) :: gl, gu
    class(graph), intent(in) :: g
    integer, intent(in) :: level
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i,j,k,n,order(2)
    integer :: Ldegrees(g%n), Udegrees(g%n)
    ! graph edge iterators
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64)

    if (level>0) then
        print *, 'ILU(k) for k>0 not implemented yet! Sorry pal.'
        call exit(1)
    endif

    order = [1, 2]
    if (present(trans)) then
        if (trans) order = [2, 1]
    endif

    allocate(cs_graph::gL)
    allocate(cs_Graph::gU)

    Ldegrees = 0
    Udegrees = 0

    ! First, determine the degrees of all the the nodes in gL, gU
    cursor = g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        call g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            if (j > i) Udegrees(i) = Udegrees(i)+1
            if (i > j) Ldegrees(i) = Ldegrees(i)+1
        enddo
    enddo

    ! Initialize gL, gU to have enough storage space for the requisite
    ! number of edges
    call gL%init(g%n,g%n,degrees=Ldegrees)
    call gU%init(g%n,g%n,degrees=Udegrees)

    ! Next, add the edges to gL, gU
    cursor = g%make_cursor(0)
    do n=1,num_blocks
        call g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            if (j>i) call gU%add_edge(i,j)
            if (i>j) call gL%add_edge(i,j)
        enddo
    enddo

    ! Finally, compress the graphs gL, gU to save on storage space and to
    ! make the immutable.
    call gL%compress()
    call gU%compress()

end subroutine incomplete_ldu_sparsity_pattern



end module ldu_solvers
