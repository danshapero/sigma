module cs_matrices

use sparse_matrices
use cs_graphs

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: cs_matrix                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    class(cs_graph), pointer :: g
    ! procedure pointers to implementations of matrix operations
    procedure(cs_find_entry_ifc), pointer, private :: find_entry
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
    procedure :: left_permute => cs_left_permute
    procedure :: right_permute => cs_right_permute
    procedure :: matvec => cs_matvec
    procedure :: matvec_t => cs_matvec_t
    procedure, private :: cs_set_value_not_preallocated
end type cs_matrix




!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    function cs_find_entry_ifc(A,i,j)
        import :: cs_matrix
        class(cs_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        integer :: cs_find_entry_ifc
    end function cs_find_entry_ifc

    subroutine cs_matvec_ifc(A,x,y)
        import :: cs_matrix, dp
        class(cs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
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
            A%find_entry => csr_find_entry

            A%matvec_impl => csr_matvec
            A%matvec_t_impl => csc_matvec

            A%left_permute_impl => cs_matrix_permute_ptrs
            A%right_permute_impl => cs_matrix_permute_vals
        case('col')
            A%find_entry => csc_find_entry

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
    integer :: k

    cs_get_value = 0_dp
    k = A%find_entry(i,j)
    if (k/=-1) cs_get_value = A%val(k)

end function cs_get_value



!--------------------------------------------------------------------------!
subroutine cs_set_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = val
    else
        call A%cs_set_value_not_preallocated(i,j,val)
    endif

end subroutine cs_set_value



!--------------------------------------------------------------------------!
subroutine cs_add_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k

    k = A%find_entry(i,j)
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%cs_set_value_not_preallocated(i,j,val)
    endif

end subroutine cs_add_value



!--------------------------------------------------------------------------!
subroutine cs_sub_matrix_add(A,B)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout)  :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k,indx

    do i=1,A%g%n
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            A%val(k) = A%val(k)+B%get_value(i,j)
        enddo
    enddo

end subroutine cs_sub_matrix_add



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

    call A%matvec_impl(x,y)

end subroutine cs_matvec



!--------------------------------------------------------------------------!
subroutine cs_matvec_t(A,x,y)                                              !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine cs_matvec_t



!--------------------------------------------------------------------------!
subroutine cs_set_value_not_preallocated(A,i,j,val)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    real(dp) :: val_temp(A%nnz)
    integer :: k

    if (A%orientation=='col') then
        call A%g%add_edge(j,i)
    else
        call A%g%add_edge(i,j)
    endif

    k = A%find_entry(i,j)

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
function csr_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: csr_find_entry

    csr_find_entry = A%g%find_edge(i,j)

end function csr_find_entry



!--------------------------------------------------------------------------!
function csc_find_entry(A,i,j)                                             !
!--------------------------------------------------------------------------!
    class(cs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: csc_find_entry

    csc_find_entry = A%g%find_edge(j,i)

end function csc_find_entry



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
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do i=1,A%g%n
        z = 0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            z = z+A%val(k)*x(j)
        enddo
        y(i) = z
    enddo

end subroutine csr_matvec



!--------------------------------------------------------------------------!
subroutine csc_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    y = 0_dp
    do j=1,A%g%n
        z = x(j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
            y(i) = y(i)+A%val(k)*z
        enddo
    enddo

end subroutine csc_matvec





end module cs_matrices
