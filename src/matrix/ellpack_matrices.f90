module ellpack_matrices

use ellpack_graphs
use sparse_matrices

implicit none



!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: ellpack_matrix                             !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:,:)
    class(ellpack_graph), pointer :: g
    ! procedure pointers to implementations of matrix operations
    procedure(ellpack_find_entry_ifc), pointer, private :: find_entry
    procedure(ellpack_permute_ifc), pointer, private :: left_permute_impl
    procedure(ellpack_permute_ifc), pointer, private :: right_permute_impl
    procedure(ellpack_matvec_ifc), pointer, private :: matvec_impl
    procedure(ellpack_matvec_ifc), pointer, private :: matvec_t_impl
contains
    ! front-ends to matrix operations
    procedure :: init => ellpack_matrix_init
    procedure :: neighbors => ellpack_matrix_neighbors
    procedure :: get_value => ellpack_get_value
    procedure :: set_value => ellpack_set_value
    procedure :: add_value => ellpack_add_value
    procedure :: sub_matrix_add => ellpack_sub_matrix_add
    procedure :: left_permute => ellpack_left_permute
    procedure :: right_permute => ellpack_right_permute
    procedure :: matvec => ellpack_matvec
    procedure :: matvec_t => ellpack_matvec_t
    procedure, private :: ellpack_set_value_not_preallocated
end type ellpack_matrix




!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine ellpack_find_entry_ifc(A,i,j,k,l)
        import :: ellpack_matrix
        class(ellpack_matrix), intent(in) :: A
        integer, intent(in)  :: i,j
        integer, intent(out) :: k,l
    end subroutine ellpack_find_entry_ifc

    subroutine ellpack_matvec_ifc(A,x,y)
        import :: ellpack_matrix, dp
        class(ellpack_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine ellpack_matvec_ifc

    subroutine ellpack_permute_ifc(A,p)
        import :: ellpack_matrix
        clasS(ellpack_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine ellpack_permute_ifc

end interface




contains



!==========================================================================!
!==========================================================================!
!== Front-ends to matrix operations                                      ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ellpack_matrix_init(A,nrow,ncol,orientation,g)                  !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

    if (present(g)) then
        select type(g)
            class is(ellpack_graph)
                A%g => g
            class default
                print *, 'Structure graph g of ellpack matrix A must be '
                print *, ' an ellpack graph. Exiting.'
                call exit(1)
        end select

        select case(orientation)
            case('row')
                A%nrow = g%n
                A%ncol = g%m
            case('col')
                A%ncol = g%n
                A%nrow = g%m
        end select
    else
        allocate(ellpack_graph::A%g)

        select case(orientation)
            case('row')
                call A%g%init(nrow,ncol)
            case('col')
                call A%g%init(ncol,nrow)
        end select
    endif

    A%nnz = A%g%ne
    allocate(A%val(A%g%max_degree,A%g%n))
    A%max_degree = A%g%max_degree

    select case(orientation)
        case('row')
            A%find_entry => ellr_find_entry

            A%matvec_impl => ellr_matvec
            A%matvec_t_impl => ellc_matvec

            A%left_permute_impl => ellpack_matrix_permute_ptrs
            A%right_permute_impl => ellpack_matrix_permute_vals
        case('col')
            A%find_entry => ellc_find_entry

            A%matvec_impl => ellc_matvec
            A%matvec_t_impl => ellr_matvec

            A%left_permute_impl => ellpack_matrix_permute_vals
            A%right_permute_impl => ellpack_matrix_permute_ptrs
    end select

end subroutine ellpack_matrix_init



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_neighbors(A,i,nbrs)                              !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine ellpack_matrix_neighbors



!--------------------------------------------------------------------------!
function ellpack_get_value(A,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: ellpack_get_value
    ! local variables
    integer :: k,l

    ellpack_get_value = 0.0_dp
    call A%find_entry(i,j,k,l)
    if (l/=-1) ellpack_get_value = A%val(l,k)

end function ellpack_get_value



!--------------------------------------------------------------------------!
subroutine ellpack_set_value(A,i,j,val)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l

    call A%find_entry(i,j,k,l)
    if (l/=-1) then
        A%val(l,k) = val
    else
        call A%ellpack_set_value_not_preallocated(i,j,val)
    endif

end subroutine ellpack_set_value



!--------------------------------------------------------------------------!
subroutine ellpack_add_value(A,i,j,val)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l

    call A%find_entry(i,j,k,l)
    if (l/=-1) then
        A%val(l,k) = A%val(l,k)+val
    else
        call A%ellpack_set_value_not_preallocated(i,j,val)
    endif

end subroutine ellpack_add_value



!--------------------------------------------------------------------------!
subroutine ellpack_sub_matrix_add(A,B)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in)     :: B
    ! local variables
    integer :: i,j,k,l,indx

    do i=1,A%g%n
        do k=1,A%max_degree
            j = A%g%node(k,i)
            if (j/=0) then
                A%val(k,i) = A%val(k,i)+B%get_value(i,j)
            endif
        enddo
    enddo

end subroutine ellpack_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine ellpack_left_permute(A,p)                                       !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine ellpack_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_right_permute(A,p)                                      !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine ellpack_right_permute



!--------------------------------------------------------------------------!
subroutine ellpack_matvec(A,x,y)                                           !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_impl(x,y)

end subroutine ellpack_matvec



!--------------------------------------------------------------------------!
subroutine ellpack_matvec_t(A,x,y)                                         !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine ellpack_matvec_t



!--------------------------------------------------------------------------!
subroutine ellpack_set_value_not_preallocated(A,i,j,val)                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l
    real(dp), allocatable :: val_temp(:,:)

    if (A%orientation=='col') then
        call A%g%add_edge(j,i)
    else
        call A%g%add_edge(i,j)
    endif

    call A%find_entry(i,j,k,l)

    if (l<=A%max_degree) then
        A%val(l,k) = val
    else
        allocate(val_temp(A%max_degree+1,A%g%n))
        val_temp(1:A%max_degree,:) = A%val
        call move_alloc(from=val_temp, to=A%val)
        A%val(A%max_degree+1,k) = val
    endif

    A%nnz = A%nnz+1
    A%max_degree = A%g%max_degree

end subroutine ellpack_set_value_not_preallocated





!==========================================================================!
!==========================================================================!
!== Implementations of matrix operations                                 ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ellr_find_entry(A,i,j,k,l)                                      !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in)  :: i,j
    integer, intent(out) :: k,l

    k = i
    l = A%g%find_edge(i,j)

end subroutine ellr_find_entry



!--------------------------------------------------------------------------!
subroutine ellc_find_entry(A,i,j,k,l)                                      !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in)  :: i,j
    integer, intent(out) :: k,l

    k = j
    l = A%g%find_edge(j,i)

end subroutine ellc_find_entry



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_permute_vals(A,p)                                !
!--------------------------------------------------------------------------!
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine ellpack_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine ellpack_matrix_permute_ptrs(A,p)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i
    real(dp) :: val(A%max_degree,A%g%n)

    do i=1,A%g%n
        val(:,p(i)) = A%val(:,i)
    enddo

    A%val = val

    call A%g%left_permute(p)

end subroutine ellpack_matrix_permute_ptrs



!--------------------------------------------------------------------------!
subroutine ellr_matvec(A,x,y)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    do i=1,A%g%n
        z = 0.0_dp
        do k=1,A%max_degree
            j = A%g%node(k,i)
            if (j/=0) z = z+A%val(k,i)*x(j)
        enddo
        y(i) = z
    enddo

end subroutine ellr_matvec



!--------------------------------------------------------------------------!
subroutine ellc_matvec(A,x,y)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z

    y = 0.0_dp
    do j=1,A%g%n
        z = x(j)
        do k=1,A%max_degree
            i = A%g%node(k,j)
            if (i/=0) y(i) = y(i)+A%val(k,j)*z
        enddo
    enddo

end subroutine ellc_matvec



end module ellpack_matrices
