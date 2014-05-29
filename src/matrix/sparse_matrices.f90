!==========================================================================!
!==========================================================================!
module sparse_matrices                                                     !
!==========================================================================!
!==========================================================================!
!==== This module contains the definition of the basic real sparse     ====!
!==== matrix class. For more optimized formats, consider using block   ====!
!==== sparse matrices or variable-block sparse matrices.               ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use linear_operator_interface
use graphs

implicit none




!--------------------------------------------------------------------------!
type, extends(linear_operator) :: sparse_matrix                            !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(graph), pointer :: g
    integer :: nnz, max_degree, order(2)
    character(len=3) :: orientation
    logical :: pos_def
    logical :: assembled
    procedure(sparse_mat_perm_ifc), pointer, private :: left_perm_impl
    procedure(sparse_mat_perm_ifc), pointer, private :: right_perm_impl
    procedure(sparse_matvec_add_ifc), pointer, private :: matvec_add_impl
    procedure(sparse_mat_get_slice_ifc), pointer, private :: get_row_impl
    procedure(sparse_mat_get_slice_ifc), pointer, private :: get_col_impl
contains
    !--------------
    ! Constructors
    !--------------
    procedure :: init => sparse_mat_init
    ! Initialize a sparse matrix given the number of rows and columns,
    ! and the desired ordering, e.g. row- or column-major ordering.

    procedure :: add_mats => sparse_mat_add_mats
    ! Initialize a sparse matrix as the sum of two other sparse matrices.

    procedure :: multiply_mats => sparse_mat_multiply_mats
    ! Initialize a sparse matrix ast he product of two sparse matrices.


    !-----------
    ! Accessors
    !-----------
    procedure :: get_value => sparse_mat_get_value
    ! Return a matrix entry

    procedure :: get_row => sparse_mat_get_row
    ! Return a row of the matrix

    procedure :: get_column => sparse_mat_get_column
    ! Return a column of the matrix


    !----------
    ! Mutators
    !----------
    procedure :: set_value => sparse_mat_set_value
    ! Set the value of a matrix entry

    procedure :: add_value => sparse_mat_add_value
    ! Add a number to a matrix entry

    procedure :: add_mat_to_self => sparse_mat_add_mat_to_self
    ! Compute the sum A <- A+B

    procedure :: zero => sparse_mat_zero
    ! Set all entries of a sparse matrix to zero

    procedure :: left_permute => sparse_mat_leftperm
    ! Permute the rows of a sparse matrix

    procedure :: right_permute => sparse_mat_rightperm
    ! Permute the columns of a sparse matrix

    procedure :: compress => sparse_mat_compress
    ! Compress the storage of a sparse matrix. This makes the underlying
    ! connectivity graph immutable


    !------------------------------
    ! Matrix-vector multiplication
    !------------------------------
    procedure :: matvec => sparse_mat_matvec
    ! Compute the product y = A*x of a matrix and a vector

    procedure :: matvec_add => sparse_mat_matvec_add
    ! Add the product A*x to the vector y


    !-------------
    ! Destructors
    !-------------
    procedure :: destroy => sparse_mat_destroy


    !--------------------------
    ! Testing, debugging & I/O
    !--------------------------
    procedure :: write_to_file => write_sparse_matrix_to_file


    !--------------------
    ! Auxiliary routines
    !--------------------
    procedure, private :: set_value_with_reallocation


    !----------
    ! Generics
    !----------
    generic :: add => add_mats, add_mat_to_self

end type sparse_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine sparse_mat_get_slice_ifc(A,nodes,vals,k)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        integer, intent(out) :: nodes(:)
        real(dp), intent(out) :: vals(:)
        integer, intent(in) :: k
    end subroutine sparse_mat_get_slice_ifc

    subroutine sparse_mat_perm_ifc(A,p)
        import :: sparse_matrix
        class(sparse_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine sparse_mat_perm_ifc

    subroutine sparse_matvec_add_ifc(A,x,y,trans)
        import :: sparse_matrix, dp
        class(sparse_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
        logical, intent(in), optional :: trans
    end subroutine sparse_matvec_add_ifc
end interface



!--------------------------------------------------------------------------!
type :: sparse_matrix_pointer                                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), pointer :: A
end type sparse_matrix_pointer




contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_init(A,nrow,ncol,orientation,g)
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(inout), optional :: g

    ! Check if the user has provided a connectivity graph
    if (present(g)) then
        ! If so, make the matrix point to it
        A%g => g

        ! Check whether the matrix is stored in row- or column-major
        ! ordering
        select case(orientation)
            case('row')
                A%nrow = g%n
                A%ncol = g%m
            case('col')
                A%ncol = g%n
                A%nrow = g%m
        end select
    else
        ! If the user has not provided a connectivity graph, initialize the
        ! graph to a default type
        allocate(ll_graph::A%g)

        ! Set the number of left- and right-nodes for the graph to be the
        ! number of rows or columns of the matrix according to whether the
        ! matrix is row- or column-oriented
        select case(orientation)
            case('row')
                call A%g%init(nrow,ncol)
            case('col')
                call A%g%init(ncol,nrow)
        end select
    endif

    ! Increment g's reference counter
    call A%g%add_reference()

    A%nnz = A%g%ne
    allocate(A%val(A%g%capacity))
    A%val = 0.0_dp
    A%max_degree = A%g%max_degree

    ! Set the order attribute for the matrix, so that the indices are
    ! reversed if we use column-major ordering
    ! I borrowed this approach from Numpy
    select case(orientation)
        case('row')
            A%order = [1, 2]

            A%get_row_impl => sparse_mat_get_slice_contiguous
            A%get_col_impl => sparse_mat_get_slice_discontiguous

            A%left_perm_impl => sparse_mat_graph_leftperm
            A%right_perm_impl => sparse_mat_graph_rightperm
        case('col')
            A%order = [2, 1]

            A%get_row_impl => sparse_mat_get_slice_discontiguous
            A%get_col_impl => sparse_mat_get_slice_contiguous

            A%left_perm_impl => sparse_mat_graph_rightperm
            A%right_perm_impl => sparse_mat_graph_leftperm
    end select

    ! Make the matrix's implementation of matvec refer to the decompressed
    ! version of the method
    A%matvec_add_impl => sparse_matvec_add_decompressed

end subroutine sparse_mat_init



!--------------------------------------------------------------------------!
subroutine sparse_mat_add_mats(A,B,C,g,orientation)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B, C
    class(graph), pointer, intent(inout), optional :: g
    character(len=3), intent(in), optional :: orientation
    ! local variables
    integer :: i,j,k,n,nv(2)
    logical :: trans1, trans2
    integer :: num_blocks, num_returned, edges(2,64)
    type(graph_edge_cursor) :: cursor
    class(graph), pointer :: ga
    character(len=3) :: ori

    !------------------------
    ! Do some error checking

    if (B%nrow/=C%nrow .or. B%ncol/=C%ncol) then
        print *, 'Dimensions for sparse matrix sum are inconsistent.'
        print *, 'Terminating.'
        call exit(1)
    endif


    !---------------------
    ! Check optional args

    ! First, check if the user has supplied a graph to use for the
    ! connectivity structure of A
    if (present(g)) then
        ! If the graph pointer isn't associated, allocate it to a CS graph
        if (.not.associated(g)) allocate(cs_graph::g)

        ! Make our local graph pointer ga point to g
        ga => g
    else
        ! If the user hasn't supplied a graph, make the local graph pointer
        ! ga a CS graph.
        allocate(cs_graph::ga)
    endif

    ! Next, check if the user has supplied an orientation for the output
    ! matrix A; if not, default to row orientation
    ori = "row"
    if (present(orientation)) ori = orientation


    !---------------------------------
    ! Build the connectivity graph ga

    ! Decide what dimension to make ga depending on the matrix orientation
    select case(ori)
        case("row")
            nv = [B%nrow, B%ncol]
            A%order = [1,2]
        case("col")
            nv = [B%ncol, B%nrow]
            A%order = [2,1]
    end select

    ! The orientation of ga is flipped relative to B%g, C%g if A has a
    ! different orientation to each matrix
    trans1 = .not.(A%order(1)==B%order(1) .and. A%order(2)==B%order(2))
    trans2 = .not.(A%order(1)==C%order(1) .and. A%order(2)==C%order(2))

    ! Compute the union of the graphs B%g, C%g and put it into ga
    call graph_union(ga,B%g,C%g,trans1,trans2)

    ! Initialize A with the connectivity structure ga
    call A%init(B%nrow,B%ncol,ori,ga)


    !-----------------------
    ! Fill the entries of A

    cursor = A%g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1
    do n=1,num_blocks
        call A%g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(A%order(1),k)
            j = edges(A%order(2),k)

            if (i/=0 .and. j/=0) then
                A%val(64*(n-1)+k) = B%get_value(i,j)+C%get_value(i,j)
            endif
        enddo
    enddo

end subroutine sparse_mat_add_mats



!--------------------------------------------------------------------------!
subroutine sparse_mat_multiply_mats(A,B,C,g,orientation)                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B, C
    class(graph), pointer, intent(inout), optional :: g
    character(len=3), intent(in), optional :: orientation
    ! local variables
    integer :: nv(2)
    class(graph), pointer :: ga, h1, h2
    character(len=3) :: ori
    logical :: tr1, tr2


    !------------------------
    ! Do some error checking

    if (B%ncol /= C%nrow) then
        print *, 'Dimensions for sparse matrix product are inconsistent.'
        print *, 'Terminating.'
        call exit(1)
    endif


    !---------------------
    ! Check optional args

    ! First, check if the user has supplied a graph to use for the
    ! connectivity structure of A
    if (present(g)) then
        ! If the graph pointer isn't associated, allocate it to a CS graph
        if (.not.associated(g)) allocate(cs_graph::g)

        ! Make our local graph pointer ga point to g
        ga => g
    else
        ! If the user hasn't supplied a graph, make the local graph pointer
        ! ga a CS graph.
        allocate(cs_graph::ga)
    endif

    ! Next, check if the user has supplied an orientation for the output
    ! matrix A; if not, default to row orientation
    ori = "row"
    if (present(orientation)) ori = orientation


    !---------------------------------
    ! Build the connectivity graph ga

    ! Decide what dimension to make ga depending on the matrix orientation
    select case(ori)
        case("row")
            nv = [B%nrow, C%ncol]
            A%order = [1,2]
        case("col")
            nv = [C%ncol, B%nrow]
            A%order = [2,1]
    end select

    ! Determine whether we're transposing matrices in the product graph
    h1 => B%g
    h2 => C%g
    tr1 = (B%orientation == "col")
    tr2 = (C%orientation == "col")

    if (ori == "col") then
        h1 => C%g
        h2 => B%g

        tr1 = xor(tr1,ori=="col")
        tr2 = xor(tr2,ori=="col")
    endif

    ! Form the graph product and put it into ga
    call graph_product(ga,h1,h2,tr1,tr2)

    ! Initialize A with the structure of ga
    call A%init(B%nrow,C%ncol,ori,ga)


    !-----------------------
    ! Fill the entries of A


end subroutine sparse_mat_multiply_mats




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function sparse_mat_get_value(A,i,j) result(val)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Set the value to return to 0
    val = 0_dp

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, return the corresponding entry in the array
    ! of values of A
    if (k/=-1) val = A%val(k)

end function sparse_mat_get_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_get_row(A,nodes,vals,k)                              !
!--------------------------------------------------------------------------!
!     This subroutine is a facade for one of either                        !
!          o sparse_mat_get_slice_contiguous                               !
!          o sparse_mat_get_slice_discontiguous.                           !
! If the matrix is stored in row-major order, then the get_row method uses !
! the implementation which assumes that row accesses are contiguous,       !
! whereas the get_column method uses the implementation assuming that      !
! column accesses require looking at every row of that column. The         !
! situation is reversed for accessing rows/columns of a matrix stored in   !
! column major order.                                                      !
!     The assignment of the function pointer get_row_impl (resp.           !
! get_col_impl) is done when initializing the matrix.                      !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: vals(:)
    integer, intent(in) :: k

    call A%get_row_impl(nodes,vals,k)

end subroutine sparse_mat_get_row



!--------------------------------------------------------------------------!
subroutine sparse_mat_get_column(A,nodes,vals,k)                           !
!--------------------------------------------------------------------------!
!     See the comment in the previous subroutine.                          !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: vals(:)
    integer, intent(in) :: k

    call A%get_col_impl(nodes,vals,k)

end subroutine sparse_mat_get_column



!--------------------------------------------------------------------------!
subroutine sparse_mat_get_slice_contiguous(A,nodes,vals,k)                 !
!--------------------------------------------------------------------------!
!     This subroutine and the one that follows it implement row/column     !
! access. This version of the algorithm assumes that the orientation of    !
! the matrix is the same as the slice being accessed, e.g. the matrix is   !
! in row-major order and we are finding all elements of a row. In that     !
! case, we can make a call to A%g%get_neighbors to find all relevant       !
! entries.                                                                 !
!     The matrix A has a function pointer which will point to either       !
! this subroutine or its discontiguous counterpart below; the function     !
! pointer is assigned at initialization of the matrix according to whether !
! the matrix is in row- or column-major order.                             !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: vals(:)
    integer, intent(in) :: k
    ! local variables
    integer :: l, d, ind

    ! Set the return values to 0
    vals = 0.0_dp

    ! Get the degree of the node k
    d = A%g%degree(k)

    ! Get all the neighbors of node k and put them in to the array `nodes`
    call A%g%get_neighbors(nodes,k)

    ! For each neighbor l of node k,
    do ind=1,d
        l = nodes(ind)

        ! get the matrix entry A(k,l)
        vals(ind) = A%val( A%g%find_edge(k,l) )
    enddo

end subroutine sparse_mat_get_slice_contiguous



!--------------------------------------------------------------------------!
subroutine sparse_mat_get_slice_discontiguous(A,nodes,vals,k)              !
!--------------------------------------------------------------------------!
!     This subroutine is the discontiguous counterpart to the subroutine   !
! above. If this subroutine is being used, it is assumed that we are       !
! accessing say an entire row of a matrix in column-major format. In that  !
! case, we cannot simply call A%g%get_neighbors in order to know which     !
! matrix entries to find -- we have to go through all entries in the row.  !
!     As you may have guessed, that is very bad and should be avoided: if  !
! the matrix you're using is in column-major ordering, you should try to   !
! only do column accesses if possible.                                     !
!     That situation can be avoided if the matrix is symmetric.            !
!     Note that for favorable access, we know that the array size needed   !
! will be at most the maximum degree of the graph. For the unfavorable     !
! access contained in this subroutine, SiGMA doesn't keep track of the     !
! maximum degree in the opposite orientation. You can easily have a graph  !
! where the maximum column degree greatly exceeds that maximum row degree  !
! or vice versa, and you could unwittingly provide arrays `nodes`, `vals`  !
! to this subroutine which have insufficient space for the operation that  !
! you're trying to do. That would result in a buffer overflow error.       !
!     TL;DR: for the love of God, try to avoid needing this subroutine.    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: vals(:)
    integer, intent(in) :: k
    ! local variables
    integer :: l, ind, next

    ! Set the index in `nodes` and `vals` of the next entry to 0
    next = 0

    ! Set the return values to 0
    vals = 0.0_dp

    ! For each node l,
    do l=1,A%g%n
        ind = A%g%find_edge(l,k)

        ! Check whether (l,k) are connected
        if (ind/=-1) then
            ! If they are,
            next = next+1
            ! put l as the next node of the slice
            nodes(next) = l
            ! and find the corresponding matrix entry.
            vals(next) = A%val(ind)
        endif
    enddo

end subroutine sparse_mat_get_slice_discontiguous



! TODO: make a version of the contiguous access subroutine which is 
! optimized for the case where A is in say CSR format, and we know that
! the entries of each row are contiguous in A%val




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_set_value(A,i,j,val)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, set the corresponding entry in the array
    ! of values of A
    if (k/=-1) then
        A%val(k) = val
    else
        call A%set_value_with_reallocation(i,j,val)
    endif

end subroutine sparse_mat_set_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_add_value(A,i,j,val)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i, j]
    ind = ind(A%order)

    ! Find the index k corresponding to the edge twixt (i,j) in A%g
    k = A%g%find_edge(ind(1),ind(2))

    ! If that edge exists, set the corresponding entry in the array
    ! of values of A
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%set_value_with_reallocation(i,j,val)
    endif

end subroutine sparse_mat_add_value



!--------------------------------------------------------------------------!
subroutine sparse_mat_add_mat_to_self(A,B)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in)    :: B
    ! local variables
    integer :: i,j,k,l,n,ind(2)
    integer :: edges(2,64), num_blocks, num_returned
    type(graph_edge_cursor) :: cursor
    type(dynamic_array) :: stack

    ! Initialize an empty stack which will store entries of B that are zero
    ! in A. For these entries, we may need to allocate additional storage
    ! space in the structure of A.
    call stack%init()

    ! Make a cursor for iterating through graph edges.
    cursor = B%g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    ! Iterate through all the non-zero entries of B.
    do n=1,num_blocks
        ! Get a chunk of edges.
        call B%g%get_edges(edges,cursor,64,num_returned)

        ! For each edge,
        do k=1,num_returned
            i = edges(B%order(1),k)
            j = edges(B%order(2),k)

            ind = [i,j]
            ind = ind(A%order)

            ! if that edge is non-null,
            if (ind(1)/=0 .and. ind(2)/=0) then
                ! find the corresponding edge in A.
                l = A%g%find_edge(ind(1),ind(2))

                ! If that edge exists in A,
                if (l/=-1) then
                    ! add the value from B to A.
                    A%val(l) = A%val(l)+B%val(64*(n-1)+k)
                else
                    ! Otherwise, B has a non-zero entry where A does not.
                    ! Push those entries onto the stack so we can expand the
                    ! storage of A later.
                    call stack%push(64*(n-1)+k)
                endif
            endif
        enddo
    enddo

    ! For now, if B is not structurally a sub-matrix of A, we throw an error
    ! because I don't feel like writing this yet.
    if (stack%length>0) then
        print *, 'Attempted to add a matrix B into a matrix A which would'
        print *, 'require expanding the storage of A. Not implemented yet.'
        call exit(1)
    endif

end subroutine sparse_mat_add_mat_to_self



!--------------------------------------------------------------------------!
subroutine sparse_mat_zero(A)                                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A

    A%val = 0_dp

end subroutine sparse_mat_zero



!--------------------------------------------------------------------------!
subroutine sparse_mat_leftperm(A,p)                                        !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_perm_impl(p)

end subroutine sparse_mat_leftperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_rightperm(A,p)                                       !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_perm_impl(p)

end subroutine sparse_mat_rightperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_graph_leftperm(A,p)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call A%g%left_permute(p,edge_p)

    ! If the edges require permutation,
    if (size(edge_p,2)>0) then
        ! Make a temporary array for the matrix values
        allocate(val(A%g%capacity))
        val = 0.0_dp

        ! Make an array of the re-arranged matrix values
        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = A%val(source:source+num_nodes-1)
        enddo

        ! Transfer the allocation of the temporary array to the array
        ! of A's matrix entries
        call move_alloc(from=val, to=A%val)
    endif

    ! Deallocate the array describing the edge permutation
    deallocate(edge_p)

end subroutine sparse_mat_graph_leftperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_graph_rightperm(A,p)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    call A%g%right_permute(p,edge_p)

    if (size(edge_p,2)>0) then
        allocate(val(A%g%capacity))
        val = 0.0_dp

        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = A%val(source:source+num_nodes-1)
        enddo

        call move_alloc(from=val, to=A%val)
    endif

    deallocate(edge_p)


end subroutine sparse_mat_graph_rightperm



!--------------------------------------------------------------------------!
subroutine sparse_mat_compress(A)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    ! local variables
    integer :: i,source,dest,num
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val(:)

    call A%g%compress(edge_p)

    if (size(edge_p,1)/=0) then
        ! Make a new array for the compressed matrix entries
        allocate(val(A%g%capacity))

        ! Move the uncompressed matrix entries into the compressed array
        do i=1,size(edge_p,1)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num = edge_p(3,i)

            val(dest:dest+num-1) = A%val(source:source+dest-1)
        enddo

        ! Transfer the allocation from the new array to the array that
        ! A points to
        call move_alloc(from=val, to=A%val)
    endif

    ! Change the strategy for the matrix's matvec to the implementation for
    ! compressed matrices, for which one needn't check
    A%matvec_add_impl => sparse_matvec_add_compressed

end subroutine sparse_mat_compress




!==========================================================================!
!==== Matrix-vector multiplication                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec(A,x,y,trans)                                  !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    logical, intent(in), optional :: trans

    y = 0_dp
    call A%matvec_add_impl(x,y,trans)

end subroutine sparse_mat_matvec



!--------------------------------------------------------------------------!
subroutine sparse_mat_matvec_add(A,x,y,trans)                              !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans

    call A%matvec_add_impl(x,y,trans)

end subroutine sparse_mat_matvec_add



!--------------------------------------------------------------------------!
subroutine sparse_matvec_add_compressed(A,x,y,trans)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i,j,k,n,order(2),edges(2,64),num_returned,num_blocks
    type(graph_edge_cursor) :: cursor

    ! Reverse the matrix orientation if we're computing A'*x instead of A*x
    order = A%order
    if (present(trans)) then
        if (trans) order = [A%order(2),A%order(1)]
    endif

    ! set the parameters for the cursor
    cursor = A%g%make_cursor(0)

    ! Set the number of blocks in which we will be performing the matrix
    ! multiplication
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        ! Get the next 64 edges of the graph
        call A%g%get_edges(edges,cursor,64,num_returned)

        ! Go through all the edges (i,j) that we just plucked from A%g
        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            ! Add A(i,j)*x(j) to y(i). Since the underlying graph is
            ! compressed, we know that i,j are not zero.
            y(i) = y(i)+A%val(64*(n-1)+k)*x(j)
        enddo
    enddo

end subroutine sparse_matvec_add_compressed



!--------------------------------------------------------------------------!
subroutine sparse_matvec_add_decompressed(A,x,y,trans)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i,j,k,n,order(2),edges(2,64),num_returned,num_blocks
    type(graph_edge_cursor) :: cursor

    ! Reverse the matrix orientation if we're computing A'*x instead of A*x
    order = A%order
    if (present(trans)) then
        if (trans) order = [A%order(2),A%order(1)]
    endif

    ! set the parameters for the cursor
    cursor = A%g%make_cursor(0)

    ! Set the number of blocks in which we will be performing the matrix
    ! multiplication
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        ! Get the next 64 edges of the graph
        call A%g%get_edges(edges,cursor,64,num_returned)

        ! Go through all the edges (i,j) that we just plucked from A%g
        do k=1,num_returned
            i = edges(order(1),k)
            j = edges(order(2),k)

            ! Check to make sure that i,j are not zero, then add 
            ! A(i,j)*x(j) to y(i); a graph which has not yet been 
            ! compressed might have null edges
            if (i/=0 .and. j/=0) then
                y(i) = y(i)+A%val(64*(n-1)+k)*x(j)
            endif
        enddo
    enddo

end subroutine sparse_matvec_add_decompressed




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_mat_destroy(A)                                           !
!--------------------------------------------------------------------------!
    class(sparse_matrix), intent(inout) :: A

    ! Deallocate the array of A's matrix entries
    deallocate(A%val)

    ! Decrement the reference counter for A%g
    call A%g%remove_reference()

    ! Nullify A's pointer to its graph. Don't de-allocate it -- there might
    ! still be other references to it someplace else.
    nullify(A%g)

end subroutine sparse_mat_destroy




!==========================================================================!
!==== Testing, debugging and I/O                                       ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine write_sparse_matrix_to_file(A,filename)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i,j,k,n
    real(dp) :: val
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64)

    open(unit=10, file=trim(filename))

    write(10,*) A%g%n, A%g%m, A%g%capacity

    cursor = A%g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1
    do n=1,num_blocks
        call A%g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            i = edges(A%order(1),k)
            j = edges(A%order(2),k)

            val = A%val(64*(n-1)+k)

            write(10,*) i,j,val
        enddo
    enddo

    close(10)

end subroutine write_sparse_matrix_to_file




!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine set_value_with_reallocation(A,i,j,val)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    class(graph), pointer :: g
    real(dp), allocatable :: vals(:)
    integer :: k, n, indx, edges(2,64), ind(2), num_blocks, num_returned
    type(graph_edge_cursor) :: cursor

    allocate(g, mold=A%g)
    call g%copy(A%g)

    ind = [i,j]
    ind = ind(A%order)

    call g%add_edge(ind(1),ind(2))

    allocate(vals(g%capacity))
    indx = g%find_edge(ind(1),ind(2))
    vals(indx) = val

    cursor = A%g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        call A%g%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            ind(1) = edges(A%order(1),k)
            ind(2) = edges(A%order(2),k)

            if (ind(1)/=0 .and. ind(2)/=0) then
                indx = g%find_edge(ind(1),ind(2))
                vals(indx) = A%val(64*(n-1)+k)
            endif
        enddo
    enddo

    call A%g%remove_reference()
    nullify(A%g)

    A%g => g
    call A%g%add_reference()
    A%nnz = A%g%ne

    call move_alloc(from=vals, to=A%val)

end subroutine set_value_with_reallocation




end module sparse_matrices
