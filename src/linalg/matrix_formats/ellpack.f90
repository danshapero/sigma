module ellpack

    use matrix
    use omp_lib

    implicit none



type, extends(sparse_matrix) :: ellpack_matrix
    integer, allocatable :: ja(:,:)
    real(kind(1d0)), allocatable :: val(:,:)
contains
    ! Constructor and accessors/mutators
    procedure :: init => ellpack_init
    procedure :: build => ellpack_build
    procedure :: get_value => ellpack_get_value
    procedure :: get_values => ellpack_get_values
    procedure :: get_neighbors => ellpack_get_neighbors
    procedure :: set_value => ellpack_set_value, & 
        & add_value => ellpack_add_value
    procedure :: set_values => ellpack_set_values, &
        & add_values => ellpack_add_values
    procedure :: zero => ellpack_zero
    procedure :: permute => ellpack_permute
    procedure :: subset_matrix_add => ellpack_subset_matrix_add
    ! matrix multiplication routines
    procedure :: matvec => ellpack_matvec
    ! forward- and back-solves for triangular systems
    procedure :: backsolve => ellpack_backsolve
    procedure :: forwardsolve => ellpack_forwardsolve
    ! routines for i/o and validation
    procedure :: convert_to_coo => ellpack_convert_to_coo
    ! auxiliary routines
    procedure :: sort_ja
end type ellpack_matrix



contains







!==========================================================================!
!==========================================================================!
!===== Accessors & mutators                                            ====!
!==========================================================================!
!==========================================================================!



!--------------------------------------------------------------------------!
subroutine ellpack_init(A,nrow,ncol,nnz,rows,cols,params)                  !
!--------------------------------------------------------------------------!
    implicit none
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:),params(:)

    A%nrow = nrow
    A%ncol = ncol
    A%nnz = nnz

    if (present(params)) then
        A%max_degree = params(1)
    else
        A%max_degree = int( ceiling(float(nnz)/nrow) )
    endif

    allocate( A%ja(A%max_degree,nrow), A%val(A%max_degree,nrow) )
    A%val = 0.d0
    A%ja = 0

    if (present(rows).and.present(cols)) then
        call A%build(rows,cols)
    endif

end subroutine ellpack_init



!--------------------------------------------------------------------------!
subroutine ellpack_build(A,rows,cols)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    ! local variables
    integer :: i,j

    do i=1,A%nnz
        do j=1,A%max_degree
            if (A%ja(j,rows(i))==0) then
                A%ja(j,rows(i)) = cols(i)
                exit
            endif
        enddo
    enddo

    call A%sort_ja()

end subroutine ellpack_build



!--------------------------------------------------------------------------!
function ellpack_get_value(A,i,j)                                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)) :: ellpack_get_value
    ! local variables
    integer :: k

    ellpack_get_value = 0.d0

    do k=1,A%max_degree
        if (A%ja(k,i)==j) ellpack_get_value = A%val(k,i)
    enddo

end function ellpack_get_value



!--------------------------------------------------------------------------!
function ellpack_get_values(A,rows,cols)                                   !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: ellpack_get_values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    ellpack_get_values = 0.d0

    ! Make this more efficient somehow
    do j=1,size(cols)
        do i=1,size(rows)
            do k=1,A%max_degree
                if (A%ja(k,rows(i))==cols(j)) then
                    ellpack_get_values(i,j) = A%val(k,i)
                endif
            enddo
        enddo
    enddo

end function ellpack_get_values



!--------------------------------------------------------------------------!
function ellpack_get_neighbors(A,row)                                      !
!--------------------------------------------------------------------------!
    implicit none
    class(ellpack_matrix), intent(in) :: A
    integer, intent(in) :: row
    integer :: ellpack_get_neighbors( A%max_degree )

    ellpack_get_neighbors = A%ja(:,row)

end function ellpack_get_neighbors



!--------------------------------------------------------------------------!
subroutine ellpack_set_value(A,i,j,val)                                    !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: k

    do k=1,A%max_degree
        if (A%ja(k,i)==j) A%val(k,i) = val
    enddo

end subroutine ellpack_set_value



!--------------------------------------------------------------------------!
subroutine ellpack_add_value(A,i,j,val)                                    !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: k

    do k=1,A%max_degree
        if (A%ja(k,i)==j) A%val(k,i) = A%val(k,i)+val
    enddo

end subroutine ellpack_add_value



!--------------------------------------------------------------------------!
subroutine ellpack_set_values(A,rows,cols,vals)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=1,A%max_degree
                if (A%ja(k,rows(i))==cols(j)) A%val(k,rows(i)) = vals(i,j)
            enddo
        enddo
    enddo

end subroutine ellpack_set_values



!--------------------------------------------------------------------------!
subroutine ellpack_add_values(A,rows,cols,vals)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k

    do j=1,size(cols)
        do i=1,size(rows)
            do k=1,A%max_degree
                if (A%ja(k,rows(i))==cols(j)) then
                    A%val(k,rows(i)) = A%val(k,rows(i))+vals(i,j)
                endif
            enddo
        enddo
    enddo

end subroutine ellpack_add_values



!--------------------------------------------------------------------------!
subroutine ellpack_zero(A)                                                 !
!--------------------------------------------------------------------------!
    implicit none
    class(ellpack_matrix), intent(inout) :: A

    A%val = 0.d0
    A%symmetric = .false.
    A%pos_def = .false.
    A%m_matrix = .false.
    A%diag_dominant = .false.

end subroutine ellpack_zero



!--------------------------------------------------------------------------!
subroutine ellpack_permute(A,p)                                            !
!--------------------------------------------------------------------------!
    implicit none
    class(ellpack_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    type(ellpack_matrix) :: B
    integer :: i,j,k,q(A%nrow)
    real(kind(1d0)) :: Aij

    allocate( B%ja(A%max_degree,A%nrow), B%val(A%max_degree,A%nrow) )
    B%ja = A%ja
    B%val = A%val
    B%max_degree = A%max_degree

    do i=1,A%nrow
        A%ja(:,p(i)) = B%ja(:,i)
        A%val(:,p(i)) = B%val(:,i)
        do k=1,A%max_degree
            if (A%ja(k,p(i))>0) A%ja(k,p(i)) = p( A%ja(k,p(i)) )
        enddo
    enddo

    ! Sort the array ja so that ja(:,i) is in order for each i
    call A%sort_ja()
   

end subroutine ellpack_permute



!--------------------------------------------------------------------------!
subroutine ellpack_subset_matrix_add(A,B)                                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    class(ellpack_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k,l

    do i=1,A%nrow
        do k=1,A%max_degree
            j = A%ja(k,i)
            if (j>0) then
                do l=1,B%max_degree
                    if (B%ja(l,i)==j) A%val(k,i) = A%val(k,i)+B%val(k,i)
                enddo
            endif
        enddo
    enddo

end subroutine ellpack_subset_matrix_add











!==========================================================================!
!==========================================================================!
!===== Matrix-vector multiplication routines                           ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine ellpack_matvec(A,x,y,rows,cols)                                 !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(out) :: y(:)
    integer, intent(in), optional :: rows(2),cols(2)
    ! local variables
    real(kind(1d0)) :: z
    integer :: i,j,k,r(2),c(2)

    r = [1,A%nrow]
    c = [1,A%ncol]
    if (present(rows)) r = rows
    if (present(cols)) c = cols

    do i=r(1),r(2)
        z = 0.d0
        do k=1,A%max_degree
            j = A%ja(k,i)
            if ( c(1)<=j .and. j<=c(2) ) z = z+A%val(k,i)*x(j)
        enddo
    enddo

end subroutine ellpack_matvec



!--------------------------------------------------------------------------!
subroutine ellpack_backsolve(A,x)                                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j,k
    real(kind(1d0)) :: Aii,z

    do i=A%nrow,1,-1
        z = x(i)
        do k=1,A%max_degree
            j = A%ja(k,i)
            if (j>i) then
                z = z-A%val(k,i)*x(j)
            elseif (j==i) then
                Aii = A%val(k,i)
            endif
        enddo
        x(i) = z/Aii
    enddo

end subroutine ellpack_backsolve



!--------------------------------------------------------------------------!
subroutine ellpack_forwardsolve(A,x)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j,k
    real(kind(1d0)) :: Aii,z

    do i=1,A%nrow
        z = x(i)
        do k=1,A%max_degree
            j = A%ja(k,i)
            if (j<i .and. j>0) then
                z = z-A%val(k,i)*x(j)
            elseif (j==i) then
                Aii = A%val(k,i)
            endif
        enddo
        x(i) = z/Aii
    enddo

end subroutine ellpack_forwardsolve









!==========================================================================!
!==========================================================================!
!======= i/o and validation routines                                   ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_convert_to_coo(A,rows,cols,vals)                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(in) :: A
    integer, intent(out) :: rows(:),cols(:)
    real(kind(1d0)), intent(out), optional :: vals(:)
    ! local variables
    integer :: i,j,k,next

    next = 0
    do i=1,A%nrow
        do k=1,A%max_degree
            j = A%ja(k,i)
            if (j/=0) then
                next = next+1
                rows(next) = i
                cols(next) = j
                if (present(vals)) vals(next) = A%val(k,i)
            endif
        enddo
    enddo

end subroutine ellpack_convert_to_coo





!==========================================================================!
!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sort_ja(A)                                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(ellpack_matrix), intent(inout) :: A
    ! local variables
    integer :: i,j,k,l,nnz_this_row
    real(kind(1d0)) :: Aij

    do i=1,A%nrow
        do k=1,A%max_degree
            if (A%ja(k,i)==0) then
                nnz_this_row = k-1
                exit
            endif
        enddo

        do k=2,nnz_this_row
            j = A%ja(k,i)
            Aij = A%val(k,i)
            do l=k-1,1,-1
                if (A%ja(l,i)<=j) exit
                A%ja(l+1,i) = A%ja(l,i)
                A%val(l+1,i) = A%val(l,i)
            enddo
            A%ja(l+1,i) = j
            A%val(l+1,i) = Aij
        enddo
    enddo

end subroutine sort_ja



end module ellpack





