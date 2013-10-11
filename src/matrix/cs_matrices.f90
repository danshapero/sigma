module cs_matrices

use sparse_matrices
use cs_graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: cs_matrix                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(cs_graph), pointer :: g
    integer :: order(2)
    ! procedure pointers to implementations of matrix operations
    procedure(cs_permute_ifc), pointer, private :: left_permute_impl
    procedure(cs_permute_ifc), pointer, private :: right_permute_impl
    procedure(cs_matvec_ifc), pointer, private :: matvec_impl
    procedure(cs_matvec_ifc), pointer, private :: matvec_t_impl
contains
    ! front-ends to matrix operations
    procedure :: init => cs_matrix_init
    procedure :: neighbors => cs_matrix_neighbors
    procedure :: get_value => cs_get_value
    procedure :: set_value => cs_set_value
    procedure :: add_value => cs_add_value
    procedure :: sub_matrix_add => cs_sub_matrix_add
    procedure :: zero => cs_zero
    procedure :: left_permute => cs_left_permute
    procedure :: right_permute => cs_right_permute
    procedure :: matvec => cs_matvec
    procedure :: matvec_t => cs_matvec_t
    procedure :: matvec_add => cs_matvec_add
    procedure :: matvec_t_add => cs_matvec_t_add
    procedure, private :: cs_set_value_not_preallocated
end type cs_matrix




!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine cs_matvec_ifc(A,x,y)
        import :: cs_matrix, dp
        class(cs_matrix), intent(in) :: A
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine

    subroutine cs_permute_ifc(A,p)
        import :: cs_matrix
        class(cs_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine cs_permute_ifc

end interface




contains



!==========================================================================!
!==========================================================================!
!== Front-ends to matrix operations                                      ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine cs_matrix_init(A,nrow,ncol,orientation,g)                       !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

    ! Check if the user has supplied a graph representing the matrix's
    ! non-zero structure
    if (present(g)) then
        ! Check to make sure the graph given is of the right type
        select type(g)
            class is(cs_graph)
                A%g => g
            class default
                ! Change this to convert the graph
                print *, 'Structure graph g of CS matrix A must be '
                print *, 'a CS graph. Exiting.'
                call exit(1)
        end select

        ! Set the number of rows and columns of the matrix to be the
        ! number of left- or right-nodes of the graph according to whether
        ! the matrix is row- or column-oriented
        select case(orientation)
            case('row')
                A%nrow = g%n
                A%ncol = g%m
            case('col')
                A%ncol = g%n
                A%nrow = g%m
        end select
    else
        ! If the user has provided no matrix structure already, allocate it
        allocate(cs_graph::A%g)

        ! Set the number of left- and right-nodes of the graph to be the
        ! number of rows or columns of the matrix according to whether the
        ! matrix is row- or column-oriented
        select case(orientation)
            case('row')
                call A%g%init(nrow,ncol)
            case('col')
                call A%g%init(ncol,nrow)
        end select
    endif

    ! Set the number of non-zero entries and the max degree of the matrix
    A%nnz = A%g%ne
    allocate(A%val(A%nnz))
    A%max_degree = A%g%max_degree

    ! Associate some procedure pointers according to whether the matrix is
    ! row- or column-oriented
    select case(orientation)
        case('row')
            A%order = [1, 2]

            A%matvec_impl => csr_matvec
            A%matvec_t_impl => csc_matvec

            A%left_permute_impl => cs_matrix_permute_ptrs
            A%right_permute_impl => cs_matrix_permute_vals
        case('col')
            A%order = [2, 1]

            A%matvec_impl => csc_matvec
            A%matvec_t_impl => csr_matvec

            A%left_permute_impl => cs_matrix_permute_vals
            A%right_permute_impl => cs_matrix_permute_ptrs
    end select

end subroutine cs_matrix_init



!--------------------------------------------------------------------------!
subroutine cs_matrix_neighbors(A,i,nbrs)                                   !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine cs_matrix_neighbors



!--------------------------------------------------------------------------!
function cs_get_value(A,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: cs_get_value
    ! local variables
    integer :: k, ind(2)

    ! (i,j) => (j,i) if the matrix is in column format
    ind = [i,j]
    ind = ind(A%order)

    ! Set the value to return to 0
    cs_get_value = 0_dp

    ! Iterate through the non-zero entries in the row/column according to
    ! the matrix's orientation
    do k=A%g%ptr(ind(1)),A%g%ptr(ind(1)+1)-1
        if (A%g%node(k)==ind(2)) cs_get_value = A%val(k)
    enddo

end function cs_get_value



!--------------------------------------------------------------------------!
subroutine cs_set_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)
    logical :: found

    ! (i,j) => (j,i) if A is in column format
    ind = [i,j]
    ind = ind(A%order)

    ! Assume that entry (i,j) is not in the matrix
    found = .false.

    ! Iterate through the non-zero entries in the row/column
    do k=A%g%ptr(ind(1)),A%g%ptr(ind(1)+1)-1
        if (A%g%node(k)==ind(2)) then
            A%val(k) = val
            found = .true.
        endif
    enddo

    ! If entry (i,j) isn't in the matrix, rebuild the matrix structure
    ! Check the indices here
    if (.not.found) call A%cs_set_value_not_preallocated(i,j,val)

end subroutine cs_set_value



!--------------------------------------------------------------------------!
subroutine cs_add_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)
    logical :: found

    ind = [i,j]
    ind = ind(A%order)

    found = .false.

    do k=A%g%ptr(ind(1)),A%g%ptr(ind(1)+1)-1
        if (A%g%node(k)==ind(2)) then
            A%val(k) = A%val(k)+val
            found = .true.
        endif
    enddo

    ! Check the indices here
    if (.not.found) call A%cs_set_value_not_preallocated(i,j,val)

end subroutine cs_add_value



!--------------------------------------------------------------------------!
subroutine cs_sub_matrix_add(A,B)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout)  :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k

    do i=1,A%g%n
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            A%val(k) = A%val(k)+B%get_value(i,j)
        enddo
    enddo

end subroutine cs_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine cs_zero(A)                                                      !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A

    A%val = 0.0_dp

end subroutine cs_zero



!--------------------------------------------------------------------------!
subroutine cs_left_permute(A,p)                                            !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine cs_left_permute



!--------------------------------------------------------------------------!
subroutine cs_right_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine cs_right_permute



!--------------------------------------------------------------------------!
subroutine cs_matvec(A,x,y)                                                !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_impl(x,y)

end subroutine cs_matvec



!--------------------------------------------------------------------------!
subroutine cs_matvec_t(A,x,y)                                              !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_t_impl(x,y)

end subroutine cs_matvec_t



!--------------------------------------------------------------------------!
subroutine cs_matvec_add(A,x,y)                                            !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_impl(x,y)

end subroutine cs_matvec_add



!--------------------------------------------------------------------------!
subroutine cs_matvec_t_add(A,x,y)                                          !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine cs_matvec_t_add



!--------------------------------------------------------------------------!
subroutine cs_set_value_not_preallocated(A,i,j,val)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    real(dp) :: val_temp(A%nnz)
    integer :: k, ind(2)

    ind = [i,j]
    ind = ind(A%order)

    call A%g%add_edge(ind(1),ind(2))

    do k=A%g%ptr(ind(1)),A%g%ptr(ind(1)+1)-1
        if (A%g%node(k)==ind(2)) exit
    enddo

    val_temp = A%val
    deallocate(A%val)
    allocate(A%val(A%nnz+1))
    A%val(1:k-1) = val_temp(1:k-1)
    A%val(k) = val
    A%val(k+1:A%nnz+1) = val_temp(k:A%nnz)

    A%nnz = A%nnz+1
    A%max_degree = A%g%max_degree

end subroutine cs_set_value_not_preallocated





!==========================================================================!
!==========================================================================!
!== Implementations of matrix operations                                 ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine cs_matrix_permute_vals(A,p)                                     !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine cs_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine cs_matrix_permute_ptrs(A,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,k,ptr(A%g%n+1)
    real(dp) :: val(A%nnz)

    do i=1,A%g%n
        ptr(p(i)+1) = A%g%ptr(i+1)-A%g%ptr(i)
    enddo

    ptr(1) = 1
    do i=1,A%g%n
        ptr(i+1) = ptr(i+1)+ptr(i)
    enddo

    do i=1,A%g%n
        do k=0,A%g%ptr(i+1)-A%g%ptr(i)-1
            val( ptr(p(i))+k ) = A%val( A%g%ptr(i)+k )
        enddo
    enddo

    A%val = val

    call A%g%left_permute(p)

end subroutine cs_matrix_permute_ptrs



!--------------------------------------------------------------------------!
subroutine csr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do i=1,A%g%n
        z = 0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            z = z+A%val(k)*x(j)
        enddo
        y(i) = y(i)+z
    enddo

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csc_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do j=1,A%g%n
        z = x(j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
            y(i) = y(i)+A%val(k)*z
        enddo
    enddo

end subroutine csc_matvec





end module cs_matrices
