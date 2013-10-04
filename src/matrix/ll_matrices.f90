module ll_matrices

use sparse_matrices
use ll_graphs
use types

implicit none


!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: ll_matrix                                  !
!--------------------------------------------------------------------------!
    real(dp), allocatable :: val(:)
    type(linked_list), allocatable :: ptrs(:)
    integer :: last, order(2)
    class(ll_graph), pointer :: g
    ! procedure pointers to implementations of matrix operations
    procedure(ll_permute_ifc), pointer, private    :: left_permute_impl
    procedure(ll_permute_ifc), pointer, private    :: right_permute_impl
    procedure(ll_matvec_ifc), pointer, private     :: matvec_impl
    procedure(ll_matvec_ifc), pointer, private     :: matvec_t_impl
contains
    ! front-ends to matrix operations
    procedure :: init => ll_matrix_init
    procedure :: neighbors => ll_matrix_neighbors
    procedure :: get_value => ll_get_value
    procedure :: set_value => ll_set_value
    procedure :: add_value => ll_add_value
    procedure :: sub_matrix_add => ll_sub_matrix_add
    procedure :: left_permute => ll_left_permute
    procedure :: right_permute => ll_right_permute
    procedure :: matvec => ll_matvec
    procedure :: matvec_t => ll_matvec_t
    procedure, private :: ll_set_value_not_preallocated
end type ll_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine ll_matvec_ifc(A,x,y)
        import :: ll_matrix, dp
        class(ll_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine

    subroutine ll_permute_ifc(A,p)
        import :: ll_matrix
        class(ll_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine ll_permute_ifc

end interface




contains



!==========================================================================!
!==========================================================================!
!== Front-ends to matrix operations                                      ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ll_matrix_init(A,nrow,ncol,orientation,g)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g
    ! local variables
    integer :: i,k

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

    if (present(g)) then
        select type(g)
            class is(ll_graph)
                A%g => g
            class default
                print *, 'Structure graph g of LL matrix A must be '
                print *, 'a LL graph. Exiting.'
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
        allocate(ll_graph::A%g)

        select case(orientation)
            case('row')
                call A%g%init(nrow,ncol)
            case('col')
                call A%g%init(ncol,nrow)
        end select
    endif

    A%nnz = A%g%ne
    allocate(A%val(A%nnz))
    A%val = 0.0_dp
    A%max_degree = A%g%max_degree

    select case(orientation)
        case('row')
            A%order = [1,2]

            A%matvec_impl   => llr_matvec
            A%matvec_t_impl => llc_matvec

            A%left_permute_impl => ll_matrix_permute_ptrs
            A%right_permute_impl => ll_matrix_permute_vals
        case('col')
            A%order = [2,1]

            A%matvec_impl   => llc_matvec
            A%matvec_t_impl => llr_matvec

            A%left_permute_impl => ll_matrix_permute_vals
            A%right_permute_impl => ll_matrix_permute_ptrs
    end select

    A%last = 1
    allocate(A%ptrs(g%n))
    do i=1,A%g%n
        do k=1,A%g%lists(i)%length
            call A%ptrs(i)%prepend(A%last)
            A%last = A%last+1
        enddo
    enddo

end subroutine ll_matrix_init



!--------------------------------------------------------------------------!
subroutine ll_matrix_neighbors(A,i,nbrs)                                   !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)

    call A%g%neighbors(i,nbrs)

end subroutine ll_matrix_neighbors



!--------------------------------------------------------------------------!
function ll_get_value(A,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: ll_get_value
    ! local variables
    integer :: k, ind(2)

    ind = [i,j]
    ind = ind(A%order)

    ll_get_value = 0_dp
    k = A%g%find_edge(ind(1),ind(2))
    if (k/=-1) ll_get_value = A%val(k)

end function ll_get_value



!--------------------------------------------------------------------------!
subroutine ll_set_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ind = [i,j]
    ind = ind(A%order)

    k = A%g%find_edge(ind(1),ind(2))
    if (k/=-1) then
        A%val(k) = val
    else
        call A%ll_set_value_not_preallocated(i,j,val)
    endif

end subroutine ll_set_value



!--------------------------------------------------------------------------!
subroutine ll_add_value(A,i,j,val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)

    ind = [i,j]
    ind = ind(A%order)

    k = A%g%find_edge(ind(1),ind(2))
    if (k/=-1) then
        A%val(k) = A%val(k)+val
    else
        call A%ll_set_value_not_preallocated(i,j,val)
    endif

end subroutine ll_add_value



!--------------------------------------------------------------------------!
subroutine ll_sub_matrix_add(A,B)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout)  :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k
    real(dp) :: Bij

    do i=1,A%g%n
        do k=1,A%g%lists(i)%length
            j = A%g%lists(i)%get_value(k)
            Bij = B%get_value(i,j)
            call A%add_value(i,j,Bij)
        enddo
    enddo

end subroutine ll_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine ll_left_permute(A,p)                                            !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine ll_left_permute



!--------------------------------------------------------------------------!
subroutine ll_right_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine ll_right_permute



!--------------------------------------------------------------------------!
subroutine ll_matvec(A,x,y)                                                !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_impl(x,y)

end subroutine ll_matvec



!--------------------------------------------------------------------------!
subroutine ll_matvec_t(A,x,y)                                              !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine ll_matvec_t



!--------------------------------------------------------------------------!
subroutine ll_set_value_not_preallocated(A,i,j,val)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k, ind(2)
    real(dp), allocatable :: val_tmp(:)

    k = A%last
    if (k>size(A%val)) then
        allocate(val_tmp(max(2*A%nnz,10)))
        val_tmp(1:A%nnz) = A%val
        call move_alloc(from=val_tmp, to=A%val)
    endif

    ind = [i,j]
    ind = ind(A%order)
    call A%g%add_edge(ind(1),ind(2))
    call A%ptrs(ind(1))%prepend(k)

    A%val(k) = val

    A%max_degree = A%g%max_degree
    A%nnz = A%nnz+1
    A%last = A%last+1

end subroutine ll_set_value_not_preallocated





!==========================================================================!
!==========================================================================!
!== Implementations of matrix operations                                 ==!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ll_matrix_permute_vals(A,p)                                     !
!--------------------------------------------------------------------------!
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine ll_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine ll_matrix_permute_ptrs(A,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,j,k
    type(linked_list) :: ptrs(A%g%n)

    do i=1,A%g%n
        do k=1,A%ptrs(i)%length
            j = A%ptrs(i)%get_value(k)
            call ptrs(i)%append(j)
        enddo
        call A%ptrs(i)%free()
    enddo

    do i=1,A%g%n
        do k=1,ptrs(i)%length
            j = ptrs(i)%get_value(k)
            call A%ptrs(p(i))%append(j)
        enddo
        call ptrs(i)%free()
    enddo

    call A%g%left_permute(p)

end subroutine ll_matrix_permute_ptrs



!--------------------------------------------------------------------------!
subroutine llr_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z,Aij

    do i=1,A%g%n
        z = 0.0_dp
        do k=1,A%ptrs(i)%length
            j = A%g%lists(i)%get_value(k)
            Aij = A%val( A%ptrs(i)%get_value(k) )
            z = z+Aij*x(j)
        enddo
        y(i) = z
    enddo

end subroutine llr_matvec



!--------------------------------------------------------------------------!
subroutine llc_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k
    real(dp) :: z,Aij

    y = 0.0_dp
    do j=1,A%g%n
        z = x(j)
        do k=1,A%ptrs(j)%length
            i = A%g%lists(j)%get_value(k)
            Aij = A%val( A%ptrs(j)%get_value(k) )
            y(i) = y(i)+Aij*z
        enddo
    enddo

end subroutine llc_matvec




end module ll_matrices
