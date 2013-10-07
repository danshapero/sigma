module block_matrices

use types
use vectors
use sparse_matrices

implicit none


!--------------------------------------------------------------------------!
type, extends(sparse_matrix) :: block_matrix                               !
!--------------------------------------------------------------------------!
    type(sparse_matrix_pointer), allocatable :: mats(:,:)
    integer :: mat_rows, mat_cols
    integer, allocatable :: row_ptr(:), col_ptr(:)

contains
    ! procedures required by virtue of inheritance from abstract
    ! sparse_matrix data type
    procedure :: init => block_mat_init
    procedure :: neighbors => block_mat_neighbors
    procedure :: get_value => block_mat_get_value
    procedure :: set_value => block_mat_set_value
    procedure :: add_value => block_mat_add_value
    procedure :: sub_matrix_add => block_mat_sub_matrix_add
    procedure :: left_permute => block_mat_left_permute
    procedure :: right_permute => block_mat_right_permute
    procedure :: matvec => block_matvec
    procedure :: matvec_t => block_matvec_t
    procedure :: matvec_add => block_matvec_add
    procedure :: matvec_t_add => block_matvec_t_add

    ! procedures specific to block matrices
    procedure :: get_block
    procedure :: set_block

    !generic :: matmul => block_mat_vector_matvec
    !generic :: matmul_t => block_mat_vector_matvec_t
end type block_matrix



contains



!--------------------------------------------------------------------------!
subroutine block_mat_init(A,nrow,ncol,orientation,g)                       !
!--------------------------------------------------------------------------!
    class(block_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

end subroutine block_mat_init



!--------------------------------------------------------------------------!
subroutine block_mat_neighbors(A,i,nbrs)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)
    ! local variables
    integer :: j

    nbrs = 0
    
end subroutine block_mat_neighbors



!--------------------------------------------------------------------------!
function block_mat_get_value(A,i,j)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: block_mat_get_value
    ! local variables
    integer :: k,l,m,n

    block_mat_get_value = 0.0_dp

    ! Find which block (k,l) contains the global entry (i,j)
    do k=1,A%mat_rows-1
        if (A%row_ptr(k)<=i .and. A%row_ptr(k+1)>i) exit
    enddo

    do l=1,A%mat_cols-1
        if (A%col_ptr(l)<=j .and. A%col_ptr(k+1)>j) exit
    enddo

    ! Local index of entry (i,j) in sub-matrix (k,l) 
    m = i-A%row_ptr(k)+1
    n = j-A%col_ptr(l)+1

    block_mat_get_value = A%mats(k,l)%A%get_value(m,n)

end function block_mat_get_value



!--------------------------------------------------------------------------!
subroutine block_mat_set_value(A,i,j,val)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l,m,n

    do k=1,A%mat_rows-1
        if (A%row_ptr(k)<=i .and. A%row_ptr(k+1)>i) exit
    enddo

    do l=1,A%mat_cols-1
        if (A%col_ptr(l)<=j .and. A%col_ptr(k+1)>j) exit
    enddo

    m = i-A%row_ptr(k)+1
    n = j-A%col_ptr(l)+1

    call A%mats(k,l)%A%set_value(m,n,val)

end subroutine block_mat_set_value



!--------------------------------------------------------------------------!
subroutine block_mat_add_value(A,i,j,val)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,l,m,n

    do k=1,A%mat_rows-1
        if (A%row_ptr(k)<=i .and. A%row_ptr(k+1)>i) exit
    enddo

    do l=1,A%mat_cols-1
        if (A%col_ptr(l)<=j .and. A%col_ptr(k+1)>j) exit
    enddo

    m = i-A%row_ptr(k)+1
    n = j-A%col_ptr(l)+1

    call A%mats(k,l)%A%add_value(m,n,val)

end subroutine block_mat_add_value



!--------------------------------------------------------------------------!
subroutine block_mat_sub_matrix_add(A,B)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables

    ! Not quite sure how I'm going to play this one

end subroutine block_mat_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine block_mat_left_permute(A,p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables

    ! This won't permute the sub-matrices themselves, it'll permute the
    ! pointers to them, changing the ordering within the block matrix

end subroutine block_mat_left_permute



!--------------------------------------------------------------------------!
subroutine block_mat_right_permute(A,p)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables

    ! Same here

end subroutine block_mat_right_permute



!--------------------------------------------------------------------------!
subroutine block_matvec(A,x,y)                                             !
!--------------------------------------------------------------------------!
    class(block_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_add(x,y)

end subroutine block_matvec



!--------------------------------------------------------------------------!
subroutine block_matvec_t(A,x,y)                                           !
!--------------------------------------------------------------------------!
    class(block_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    y = 0.0_dp
    call A%matvec_t_add(x,y)

end subroutine block_matvec_t



!--------------------------------------------------------------------------!
subroutine block_matvec_add(A,x,y)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i1,j1,i2,j2,k,l

    do l=1,A%mat_cols
        j1 = A%col_ptr(l)
        j2 = A%col_ptr(l+1)-1

        do k=1,A%mat_rows
            i1 = A%row_ptr(k)
            i2 = A%row_ptr(k+1)-1

            call A%mats(k,l)%A%matvec_add(x(j1:j2),y(i1:i2))
        enddo
    enddo

end subroutine block_matvec_add



!--------------------------------------------------------------------------!
subroutine block_matvec_t_add(A,x,y)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(in) :: A
    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: y(:)
    ! local variables
    integer :: i1,j1,i2,j2,k,l

    do l=1,A%mat_cols
        j1 = A%col_ptr(l)
        j2 = A%col_ptr(l+1)-1

        do k=1,A%mat_rows
            i1 = A%row_ptr(k)
            i2 = A%row_ptr(k+1)-1

            call A%mats(k,l)%A%matvec_t_add(x(i1:i2),y(j1:j2))
        enddo
    enddo

end subroutine block_matvec_t_add



!--------------------------------------------------------------------------!
function get_block(A,k,l)                                                  !
!--------------------------------------------------------------------------!
    class(block_matrix), intent(in) :: A
    integer, intent(in) :: k,l
    class(sparse_matrix), pointer :: get_block

    get_block => A%mats(k,l)%A

end function get_block



!--------------------------------------------------------------------------!
subroutine set_block(A,B,k,l)                                              !
!--------------------------------------------------------------------------!
    class(block_matrix), intent(inout) :: A
    class(sparse_matrix), pointer, intent(in) :: B
    integer, intent(in) :: k,l

    A%mats(k,l)%A => B

end subroutine set_block



!--------------------------------------------------------------------------!
subroutine block_mat_vector_matvec(A,x,y)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(block_matrix), intent(in) :: A
    class(vector), intent(in)  :: x
    class(vector), intent(out) :: y
    ! local variables

end subroutine block_mat_vector_matvec



end module block_matrices
