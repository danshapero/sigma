module bsr_matrix_mod

    use sparse_matrix_mod
    use omp_lib

    implicit none



type, extends(sparse_matrix) :: bsr_matrix
    integer :: mrow,mcol,nbr,nbc,nnzb
    integer, allocatable :: ia(:), ja(:)
    real(kind(1d0)), allocatable :: val(:,:,:)
contains
    ! Constructor and accessors/mutators
    procedure :: init => bsr_init
    procedure :: build => bsr_build
    procedure :: get_value => bsr_get_value
    procedure :: get_values => bsr_get_values
    procedure :: get_neighbors => bsr_get_neighbors
    procedure :: set_value => bsr_set_value, add_value => bsr_add_value
    procedure :: set_values => bsr_set_values, add_values => bsr_add_values
    procedure :: permute => bsr_permute
    procedure :: subset_matrix_add => bsr_subset_matrix_add
    ! matrix multiplication routines
    procedure :: matvec => bsr_matvec
    ! forward- and back-solves for triangular systems
    procedure :: backsolve => bsr_backsolve
    procedure :: forwardsolve => bsr_forwardsolve
    ! routines for i/o and validation
    procedure :: convert_to_coo => bsr_convert_to_coo
    procedure :: write_to_file => bsr_write_to_file
    ! auxiliary routines
    procedure :: sort_ja
    procedure :: whole_block
end type bsr_matrix



contains






!==========================================================================!
!==========================================================================!
!==== Accessors & mutators                                             ====!
!==========================================================================!
!==========================================================================!



!--------------------------------------------------------------------------!
subroutine bsr_init(A,nrow,ncol,nnz,rows,cols,vals,params)                 !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: nrow,ncol,nnz
    integer, intent(in), optional :: rows(:),cols(:),params(:)
    real(kind(1d0)), intent(in), optional :: vals(:)

    A%nrow = nrow
    A%ncol = ncol
    A%nnz = nnz

    if (present(params)) then
        A%mrow = params(1)
        if ( mod(A%mrow,nrow)/=0 .and. A%mrow/=1 ) then
            print *, 'Error in bsr_init:'
            print *, '   The row block size you have specified for a BSR'
            print *, '   matrix does not divide the number of rows of   '
            print *, '   the  matrix.                                   '
            print *, '   Reverting to a row block size of 1, for which  '
            print *, '   the BSR format is no more efficient than CSR.  '
            A%mrow = 1
        endif

        if (size(params)>1) then
            A%mcol = params(2)
            if ( mod(A%mcol,ncol)/=0 .and. A%mcol/=1 ) then
                print *, 'Error in bsr_init:'
                print *, '   The column block size specified for a BSR  '
                print *, '   matrix does not divide the number of       '
                print *, '   columns of the matrix.                     '
                print *, '   Reverting to a column block size of 1, for '
                print *, '   BSR format is no more efficient than CSR.  '
            endif
        else
            A%mcol = A%mrow
        endif
    else
        A%mrow = 1
        A%mcol = 1
    endif

    A%nbr = nrow/A%mrow
    A%nbc = ncol/A%mcol
    A%nnzb = nnz/(A%mrow * A%mcol)

    allocate( A%ia(A%mrow+1), A%ja(A%nnzb), A%val(A%mrow,A%mcol,A%nnzb) )
    A%ia = 0
    A%ia(A%mrow+1) = A%nnzb+1
    A%ja = 0

    if (present(rows).and.present(cols)) then
        call A%build(rows,cols,vals)
    endif

end subroutine bsr_init



!--------------------------------------------------------------------------!
subroutine bsr_build(A,rows,cols,vals)                                     !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in), optional :: vals(:)
    ! local variables
    integer :: i,j,k,n,startptr
    integer :: work(A%mrow)

    ! Fill in a work array: work(i) = #of non-zero blocks in (block)-row i
    work = 0
    if ( size(rows)==A%nnzb ) then
        ! If the structure provided by the arrays rows,cols corresponds to
        ! the connections between blocks and not between the raw unknowns,
        ! we can fill in the work array directly
        do i=1,A%nnzb
            work( rows(i) ) = work( rows(i) )+1
        enddo
    else
        ! If the structure provided by the arrays rows,cols coresponds to
        ! the raw unknowns and not to the blocks, we have to reduce it to
        ! the connections between the blocks

        ! Haven't done this yet. Probably need a linked list data structure
        ! to do it right.
    endif

    A%max_degree = maxval(work)

    ! Fill in the array ia; ia(i+1)-ia(i) = #of non-zero entries in row i
    A%ia(1) = 1
    do i=1,A%mrow
        A%ia(i+1) = A%ia(i)+work(i)
    enddo

    ! Fill in the array ja; ja(j) = the (block)-column of non-zero block #j
    work = 0
    if ( size(cols)==A%nnzb) then
        do i=1,A%nnz
            startptr = A%ia( rows(i) )
            A%ja( startptr+work(rows(i)) ) = cols(i)
            work( rows(i) ) = work( rows(i) )+1
        enddo
    else

    endif

    ! Sort the array ja: ja(ia(i)),ja(ia(i+1)-1) are in order
    call A%sort_ja()

end subroutine bsr_build



!--------------------------------------------------------------------------!
function bsr_get_value(A,i,j)                                              !
!--------------------------------------------------------------------------!
    implicit none
    class(bsr_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)) :: bsr_get_value
    ! local variables
    integer :: k,row,col,rb,cb

    bsr_get_value = 0.d0

    row = (i-1)/A%nbr+1
    col = (j-1)/A%nbc+1
    rb = mod(i-1,A%nbr)+1
    cb = mod(j-1,A%nbc)+1

    do k=A%ia( row ),A%ia( row+1 )-1
        if (A%ja(k) == col) then
            bsr_get_value = A%val(rb,cb,k)
        endif
    enddo

end function bsr_get_value



!--------------------------------------------------------------------------!
function bsr_get_values(A,rows,cols)                                       !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)) :: bsr_get_values(size(rows),size(cols))
    ! local variables
    integer :: i,j,k,rb,cb,row,col

    bsr_get_values = 0.d0

    if (A%whole_block(rows,cols)) then
        row = (rows(1)-1)/A%nbr+1
        col = (cols(1)-1)/A%nbc+1
        do k=A%ia(row),A%ia(row+1)-1
            if (A%ja(k)==col) then
                bsr_get_values = A%val(:,:,k)
            endif
        enddo
    else
        do j=1,size(cols)
            col = (cols(j)-1)/A%nbc+1
            cb = mod(cols(j)-1,A%nbc)+1
            do i=1,size(rows)
                row = (rows(i)-1)/A%nbr+1
                rb = mod(rows(i)-1,A%nbc)+1
                do k=A%ia(row),A%ia(row+1)-1
                    if (A%ja(k)==col) then
                        bsr_get_values(i,j) = A%val(rb,cb,k)
                    endif
                enddo
            enddo
        enddo
    endif

end function bsr_get_values



!--------------------------------------------------------------------------!
function bsr_get_neighbors(A,row)                                          !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    integer, intent(in) :: row
    integer :: bsr_get_neighbors( A%max_degree )
    ! local variables
    integer :: start,finish

    bsr_get_neighbors = 0
    start = A%ia(row)
    finish = A%ia(row+1)-1
    bsr_get_neighbors(1:finish-start+1) = A%ja(start:finish)

end function bsr_get_neighbors



!--------------------------------------------------------------------------!
subroutine bsr_set_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: row,col,rb,cb,k

    row = (i-1)/A%nbr+1
    col = (j-1)/A%nbc+1
    rb = mod(i-1,A%nbr)+1
    cb = mod(j-1,A%nbc)+1

    do k=A%ia(row),A%ia(row+1)-1
        if (A%ja(k)==col) then
            A%val(rb,cb,k) = val
        endif
    enddo

end subroutine bsr_set_value



!--------------------------------------------------------------------------!
subroutine bsr_add_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(kind(1d0)), intent(in) :: val
    ! local variables
    integer :: k,row,col,rb,cb

    row = (i-1)/A%nbr+1
    col = (j-1)/A%nbc+1
    rb = mod(i-1,A%nbr)+1
    cb = mod(j-1,A%nbc)+1

    do k=A%ia(row),A%ia(row+1)-1
        if (A%ja(k)==col) then
            A%val(rb,cb,k) = A%val(rb,cb,k)+val
        endif
    enddo

end subroutine bsr_add_value



!--------------------------------------------------------------------------!
subroutine bsr_set_values(A,rows,cols,vals)                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k,row,col,rb,cb

    if (A%whole_block(rows,cols)) then
        row = (rows(1)-1)/A%nbr+1
        col = (cols(1)-1)/A%nbc+1
        do k=A%ia(row),A%ia(row+1)-1
            if (A%ja(k)==col) A%val(:,:,k) = vals
        enddo
    else
        do j=1,size(cols)
            col = (cols(j)-1)/A%nbc+1
            cb = mod(cols(j)-1,A%nbc)+1
            do i=1,size(rows)
                row = (rows(i)-1)/A%nbr+1
                rb = mod(rows(i)-1,A%nbc)+1
                do k=A%ia(row),A%ia(row+1)-1
                    if (A%ja(k)==col) then
                        A%val(rb,cb,k) = vals(i,j)
                    endif
                enddo
            enddo
        enddo
    endif

end subroutine bsr_set_values



!--------------------------------------------------------------------------!
subroutine bsr_add_values(A,rows,cols,vals)                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: rows(:),cols(:)
    real(kind(1d0)), intent(in) :: vals(size(rows),size(cols))
    ! local variables
    integer :: i,j,k,row,col,rb,cb

    if (A%whole_block(rows,cols)) then
        row = (rows(1)-1)/A%nbr+1
        col = (cols(1)-1)/A%nbc+1
        do k=A%ia(row),A%ia(row+1)-1
            if (A%ja(k)==col) A%val(:,:,k) = A%val(:,:,k)+vals
        enddo
    else
        do j=1,size(cols)
            col = (cols(j)-1)/A%nbc+1
            cb = mod(cols(j)-1,A%nbc)+1
            do i=1,size(rows)
                row = (rows(i)-1)/A%nbr+1
                rb = mod(rows(i)-1,A%nbc)+1
                do k=A%ia(row),A%ia(row+1)-1
                    if (A%ja(k)==col) then
                        A%val(rb,cb,k) = A%val(rb,cb,k)+vals(i,j)
                    endif
                enddo
            enddo
        enddo
    endif

end subroutine bsr_add_values



!--------------------------------------------------------------------------!
subroutine bsr_permute(A,p)                                                !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    type(bsr_matrix) :: B
    integer :: i,j,k,n,astart,bstart,q( A%mrow )
    real(kind(1d0)) :: Aij(A%nbr,A%nbc)

    do i=1,A%mrow
        q( p(i) ) = i
    enddo

    allocate( B%ia( A%mrow+1 ), B%ja( A%nnzb ), &
         & B%val( A%nbr, A%nbc, A%nnzb ) )

    B%ia = A%ia
    B%ja = A%ja
    B%val = A%val

    A%ia = 0
    A%ja = 0
    A%val = 0.d0

    A%ia(1) = 1
    do i=1,A%mrow
        astart = A%ia(i)
        bstart = B%ia(q(i))
        A%ia(i+1) = astart+B%ia(q(i)+1)-B%ia(q(i))
        do j=0,A%ia(i+1)-A%ia(i)-1
            A%ja(astart+j) = p( B%ja(bstart+j) )
            A%val(:,:,astart+j) = B%val(:,:,bstart+j)
        enddo
    enddo

    ! Sort the array ja so that ja(ia(i)),ja(ia(i+1)-1) are in order
    call A%sort_ja()

end subroutine bsr_permute



!--------------------------------------------------------------------------!
subroutine bsr_subset_matrix_add(A,B)                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    class(bsr_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k

    !$omp parallel do private(j,k)
    do i=1,A%mrow
        do j=A%ia(i),A%ia(i+1)-1
            do k=B%ia(i),B%ia(i+1)-1
                if ( B%ja(k)==A%ja(j) ) then
                    A%val(:,:,j) = A%val(:,:,j)+B%val(:,:,k)
                endif
            enddo
        enddo
    enddo
    !$omp end parallel do

    A%symmetric = A%symmetric .and. B%symmetric
    A%pos_def = A%pos_def .and. B%pos_def
    A%diag_dominant = A%diag_dominant .and. B%diag_dominant

end subroutine bsr_subset_matrix_add











!==========================================================================!
!==========================================================================!
!===== Matrix-vector multiplication routines                           ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine bsr_matvec(A,x,y,rows,cols)                                     !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    real(kind(1d0)), intent(in) :: x(:)
    real(kind(1d0)), intent(out) :: y(:)
    integer, intent(in), optional :: rows(2),cols(2)
    ! local variables
    real(kind(1d0)) :: z(A%nbr)
    integer :: i,j,k,r(2),c(2)

    r = [1,A%mrow]
    c = [1,A%mcol]
    if (present(rows)) r = rows
    if (present(cols)) c = cols

    !$omp parallel do private(j,k,z)
    do i=r(1),r(2)
        z = 0.d0
        do k=A%ia(i),A%ia(i+1)-1
            j = A%ja(k)
            if ( c(1)<=j .and. j<=c(2) ) then
                z = z+matmul(A%val(:,:,k),x(A%nbc*(j-1)+1:A%nbc*j))
            endif
        enddo
        y(A%nbr*(i-1)+1:A%nbr*i) = z
    enddo
    !$omp end parallel do

end subroutine bsr_matvec



!--------------------------------------------------------------------------!
subroutine bsr_backsolve(A,x)                                              !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j,k,piv(A%nbr),info
    real(kind(1d0)) :: Aii(A%nbr,A%nbc),z(A%nbr)

    associate( nbr=>A%nbr, nbc=>A%nbc, mrow=>A%mrow, mcol=>A%mcol )

    do i=mrow,1,-1
        z = x( nbr*(i-1)+1:nbr*i )
        do j=A%ia(i),A%ia(i+1)-1
            k = A%ja(j)
            if (k>i) then
                z = z-matmul( A%val(:,:,j), x(nbc*(k-1)+1:nbc*k) )
            elseif (k==i) then
                Aii = A%val(:,:,j)
            endif
        enddo
        ! Make a call to LAPACK here
        call dgesv( nbr,1,Aii,nbr,piv,z,nbr,info)

        ! Should we prepare the Cholesky/LU factorization of the diagonal
        ! blocks in advance?
        ! It will be inefficient to repeatedly factor the blocks if we have
        ! to use a block Jacobi/Gauss-Seidel preconditioner
        ! However, we can make a variable "factorization_is_up_to_date"
        ! which gets set to true once a factorization is computed, then set
        ! to false whenever the user calls add_value/set_value
        x( nbr*(i-1)+1:nbr*i ) = z
    enddo

    end associate

end subroutine bsr_backsolve



!--------------------------------------------------------------------------!
subroutine bsr_forwardsolve(A,x)                                           !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    real(kind(1d0)), intent(inout) :: x(:)
    ! local variables
    integer :: i,j,k,piv(A%nbr),info
    real(kind(1d0)) :: Aii(A%nbr,A%nbc),z(A%nbr)

    associate( nbr=>A%nbr, nbc=>A%nbc, mrow=>A%mrow, mcol=>A%mcol )

    do i=1,mrow
        z = x( nbr*(i-1)+1:nbr*i )
        do j=A%ia(i),A%ia(i+1)-1
            k = A%ja(j)
            if (k>i) then
                z = z-matmul( A%val(:,:,j), x(nbc*(k-1)+1:nbc*k) )
            elseif (k==i) then
                Aii = A%val(:,:,j)
            endif
        enddo
        ! Make a call to LAPACK here
        call dgesv( nbr,1,Aii,nbr,piv,z,nbr,info)

        ! Should we prepare the Cholesky/LU factorization of the diagonal
        ! blocks in advance?
        ! It will be inefficient to repeatedly factor the blocks if we have
        ! to use a block Jacobi/Gauss-Seidel preconditioner
        ! However, we can make a variable "factorization_is_up_to_date"
        ! which gets set to true once a factorization is computed, then set
        ! to false whenever the user calls add_value/set_value
        x( nbr*(i-1)+1:nbr*i ) = z
    enddo

    end associate

end subroutine bsr_forwardsolve












!==========================================================================!
!==========================================================================!
!======= i/o and validation routines                                   ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine bsr_convert_to_coo(A,rows,cols,vals)                            !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    integer, intent(out) :: rows(:),cols(:)
    real(kind(1d0)), intent(out), optional :: vals(:)
    ! local variables
    integer :: i,j

end subroutine bsr_convert_to_coo


!--------------------------------------------------------------------------!
subroutine bsr_write_to_file(A,filename)                                   !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i,j

    open(unit=100,file=trim(filename)//".txt")

    close(100)

end subroutine bsr_write_to_file








!==========================================================================!
!==========================================================================!
!==== Auxiliary subroutines                                             ===!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sort_ja(A)                                                      !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(inout) :: A
    ! local variables
    integer :: i,j,k,n
    real(kind(1d0)) :: Aij(A%nbr,A%nbc)

    do i=1,A%mrow
        do j=A%ia(i)+1,A%ia(i+1)-1
            n = A%ja(j)
            Aij = A%val(:,:,j)
            do k=j-1,A%ia(i),-1
                if (A%ja(k)<=n) exit
                A%ja(k+1) = A%ja(k)
                A%val(:,:,k+1) = A%val(:,:,k)
            enddo
            A%ja(k+1) = n
            A%val(:,:,k+1) = Aij
        enddo
    enddo

end subroutine sort_ja



!--------------------------------------------------------------------------!
logical function whole_block(A,rows,cols)                                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(bsr_matrix), intent(in) :: A
    integer, intent(in) :: rows(:),cols(:)
    ! local variables
    integer :: i

    whole_block = size(rows)==A%nbr .and. size(cols)==A%nbc
    do i=1,size(rows)
        whole_block = whole_block .and. mod(rows(i)-1,A%nbr)==i-1
    enddo
    do i=1,size(cols)
        whole_block = whole_block .and. mod(cols(i)-1,A%nbc)==i-1
    enddo

end function whole_block


end module bsr_matrix_mod
