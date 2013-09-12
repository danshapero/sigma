module bcs_matrices

use sparse_matrices
use block_sparse_matrices
use cs_graphs

implicit none


!--------------------------------------------------------------------------!
type, extends(block_sparse_matrix) :: bcs_matrix                           !
!--------------------------------------------------------------------------!
    class(cs_graph), pointer :: g
    real(dp), allocatable :: val(:,:,:)
    ! procedure pointers to implementations of matrix operations
    procedure(bcs_find_block_ifc), pointer, private :: find_block
    procedure(bcs_permute_ifc), pointer, private :: left_permute_impl
    procedure(bcs_permute_ifc), pointer, private :: right_permute_impl
    procedure(bcs_matvec_ifc), pointer, private :: matvec_impl
    procedure(bcs_matvec_ifc), pointer, private :: matvec_t_impl
    procedure(bcs_block_matvec_ifc), pointer, private :: block_matvec_impl
    procedure(bcs_block_matvec_ifc), pointer, private :: block_matvec_t_impl
    procedure(bcs_l_block_matvec_ifc), pointer, private :: l_block_matvec_impl
    procedure(bcs_l_block_matvec_ifc), pointer, private :: l_block_matvec_t_impl
    procedure(bcs_r_block_matvec_ifc), pointer, private :: r_block_matvec_impl
    procedure(bcs_r_block_matvec_ifc), pointer, private :: r_block_matvec_t_impl
contains
    ! front-ends to matrix operations
    procedure :: init => bcs_init
    procedure :: neighbors => bcs_matrix_neighbors
    procedure :: get_value => bcs_get_value
    procedure :: set_value => bcs_set_value
    procedure :: add_value => bcs_add_value
    procedure :: get_block => bcs_get_block
    procedure :: set_block => bcs_set_block
    procedure :: add_block => bcs_add_block
    procedure :: sub_matrix_add => bcs_sub_matrix_add
    procedure :: left_permute => bcs_left_permute
    procedure :: right_permute => bcs_right_permute
    procedure :: matvec => bcs_matvec
    procedure :: matvec_t => bcs_matvec_t
    procedure :: block_matvec => bcs_block_matvec
    procedure :: block_matvec_t => bcs_block_matvec_t
    procedure :: l_block_matvec => bcs_l_block_matvec
    procedure :: l_block_matvec_t => bcs_l_block_matvec_t
    procedure :: r_block_matvec => bcs_r_block_matvec
    procedure :: r_block_matvec_t => bcs_r_block_matvec_t
    procedure, private :: bcs_set_block_not_preallocated
end type bcs_matrix



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    function bcs_find_block_ifc(A,i,j)
        import :: bcs_matrix
        class(bcs_matrix), intent(in) :: A
        integer, intent(in) :: i,j
        integer :: bcs_find_block_ifc
    end function bcs_find_block_ifc

    subroutine bcs_permute_ifc(A,p)
        import :: bcs_matrix
        class(bcs_matrix), intent(inout) :: A
        integer, intent(in) :: p(:)
    end subroutine bcs_permute_ifc

    subroutine bcs_matvec_ifc(A,x,y)
        import :: bcs_matrix, dp
        class(bcs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
    end subroutine bcs_matvec_ifc

    subroutine bcs_block_matvec_ifc(A,x,y)
        import :: bcs_matrix, dp
        class(bcs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:,:)
        real(dp), intent(out) :: y(:,:)
    end subroutine bcs_block_matvec_ifc

    subroutine bcs_l_block_matvec_ifc(A,x,y)
        import :: bcs_matrix, dp
        class(bcs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:,:)
    end subroutine bcs_l_block_matvec_ifc

    subroutine bcs_r_block_matvec_ifc(A,x,y)
        import :: bcs_matrix, dp
        class(bcs_matrix), intent(in) :: A
        real(dp), intent(in)  :: x(:,:)
        real(dp), intent(out) :: y(:)
    end subroutine bcs_r_block_matvec_ifc

end interface





contains



!==========================================================================!
!==========================================================================!
!==== Front-ends to matrix operations                                  ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine bcs_init(A,nrow,ncol,orientation,g)                             !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: nrow, ncol
    character(len=3), intent(in) :: orientation
    class(graph), pointer, intent(in), optional :: g

    A%nrow = nrow
    A%ncol = ncol
    A%orientation = orientation

    if (present(g)) then
        ! Check to make sure the graph given is of the right type
        select type(g)
            class is(cs_graph)
                A%g => g
            class default
                ! Change this to converting the graph
                print *, 'Structure graph g of a Block CS matrix A must be'
                print *, 'a CS graph. Exiting.'
                call exit(1)
        end select

        select case(orientation)
            case('row')
                A%nr = nrow/g%n
                A%nc = ncol/g%m
            case('col')
                A%nc = ncol/g%n
                A%nr = nrow/g%m
        end select
    else
        ! What the hell do we even do here?
    endif

    A%nnz = A%g%ne
    A%max_degree = A%g%max_degree

    select case(orientation)
        case('row')
            allocate(A%val(A%nr,A%nc,A%nnz))

            A%find_block => bcsr_find_block

            A%left_permute_impl => bcs_matrix_permute_ptrs
            A%right_permute_impl => bcs_matrix_permute_vals

            A%matvec_impl   => bcsr_matvec
            A%matvec_t_impl => bcsc_matvec

            A%block_matvec_impl   => bcsr_block_matvec
            A%block_matvec_t_impl => bcsc_block_matvec

            A%l_block_matvec_impl   => bcsr_l_block_matvec
            A%l_block_matvec_t_impl => bcsc_l_block_matvec

            A%r_block_matvec_impl   => bcsr_r_block_matvec
            A%r_block_matvec_t_impl => bcsc_r_block_matvec
        case('col')
            allocate(A%val(A%nc,A%nr,A%nnz))

            A%find_block => bcsc_find_block

            A%left_permute_impl => bcs_matrix_permute_vals
            A%right_permute_impl => bcs_matrix_permute_ptrs

            A%matvec_impl   => bcsc_matvec
            A%matvec_t_impl => bcsr_matvec

            A%block_matvec_impl   => bcsc_block_matvec
            A%block_matvec_t_impl => bcsr_block_matvec

            A%l_block_matvec_impl   => bcsc_l_block_matvec
            A%l_block_matvec_t_impl => bcsr_l_block_matvec

            A%r_block_matvec_impl   => bcsc_r_block_matvec
            A%r_block_matvec_t_impl => bcsr_r_block_matvec
    end select

end subroutine bcs_init



!--------------------------------------------------------------------------!
subroutine bcs_matrix_neighbors(A,i,nbrs)                                  !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    integer, intent(in)  :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    call A%g%neighbors(i,nbrs)

end subroutine bcs_matrix_neighbors



!--------------------------------------------------------------------------!
function bcs_get_value(A,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp) :: bcs_get_value
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    bcs_get_value = 0_dp

    k = A%find_block(rb,cb)
    if (k/=-1) then
        if (A%orientation=="col") then
            bcs_get_value = A%val(jb,ib,k)
        else
            bcs_get_value = A%val(ib,jb,k)
        endif
    endif

end function bcs_get_value



!--------------------------------------------------------------------------!
subroutine bcs_set_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    k = A%find_block(rb,cb)
    if (k/=-1) then
        if (A%orientation=="col") then
            A%val(jb,ib,k) = val
        else
            A%val(ib,jb,k) = val
        endif
    endif
    ! Change this to an else statement and add soemthing to add values in
    ! space which hasn't been preallocated yet

end subroutine bcs_set_value



!--------------------------------------------------------------------------!
subroutine bcs_add_value(A,i,j,val)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: val
    ! local variables
    integer :: k,rb,cb,ib,jb

    rb = (i-1)/A%nr+1
    cb = (j-1)/A%nc+1
    ib = i-A%nr*rb
    jb = j-A%nc*cb

    k = A%find_block(rb,cb)
    if (k/=-1) then
        if (A%orientation=="col") then
            A%val(jb,ib,k) = A%val(jb,ib,k)+val
        else
            A%val(ib,jb,k) = A%val(ib,jb,k)+val
        endif
    endif
    ! Change this to an else statement and add soemthing to add values in
    ! space which hasn't been preallocated yet

end subroutine bcs_add_value



!--------------------------------------------------------------------------!
subroutine bcs_get_block(A,i,j,vals)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    real(dp), intent(out) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%find_block(i,j)
    if (A%orientation=="col") then
        vals = transpose(A%val(:,:,k))
    else
        vals = A%val(:,:,k)
    endif

end subroutine bcs_get_block



!--------------------------------------------------------------------------!
subroutine bcs_set_block(A,i,j,vals)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%find_block(i,j)
    if (A%orientation=="col") then
        A%val(:,:,k) = transpose(vals)
    else
        A%val(:,:,k) = vals
    endif

end subroutine bcs_set_block



!--------------------------------------------------------------------------!
subroutine bcs_add_block(A,i,j,vals)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
    real(dp), intent(in) :: vals(:,:)
    ! local variables
    integer :: k

    k = A%find_block(i,j)
    if (A%orientation=="col") then
        A%val(:,:,k) = A%val(:,:,k)+transpose(vals)
    else
        A%val(:,:,k) = A%val(:,:,k)+vals
    endif

end subroutine bcs_add_block



!--------------------------------------------------------------------------!
subroutine bcs_sub_matrix_add(A,B)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    class(sparse_matrix), intent(in) :: B
    ! local variables
    integer :: i,j,k,n1,n2
    real(dp) :: vals(size(A%val,1),size(A%val,2))

    select type(B)
        class is(block_sparse_matrix)
            do i=1,A%g%n
                do k=A%g%ptr(i),A%g%ptr(i+1)-1
                    j = A%g%node(k)
                    call B%get_block(i,j,vals)
                    A%val(:,:,k) = A%val(:,:,k)+vals
                enddo
            enddo
            return
        class is(sparse_matrix)
!            do i=1,A%g%n
!                do k=A%g%ptr(i),A%g%ptr(i+1)-1
!                    j = A%g%node(k)
!                    
!                enddo
!            enddo
    end select

end subroutine bcs_sub_matrix_add



!--------------------------------------------------------------------------!
subroutine bcs_left_permute(A,p)                                           !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%left_permute_impl(p)

end subroutine bcs_left_permute



!--------------------------------------------------------------------------!
subroutine bcs_right_permute(A,p)                                          !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%right_permute_impl(p)

end subroutine bcs_right_permute



!--------------------------------------------------------------------------!
subroutine bcs_matvec(A,x,y)                                               !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_impl(x,y)

end subroutine bcs_matvec



!--------------------------------------------------------------------------!
subroutine bcs_matvec_t(A,x,y)                                             !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)

    call A%matvec_t_impl(x,y)

end subroutine bcs_matvec_t



!--------------------------------------------------------------------------!
subroutine bcs_block_matvec(A,x,y)                                         !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)

    call A%block_matvec_impl(x,y)

end subroutine bcs_block_matvec



!--------------------------------------------------------------------------!
subroutine bcs_block_matvec_t(A,x,y)                                       !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)

    call A%block_matvec_t_impl(x,y)

end subroutine bcs_block_matvec_t



!--------------------------------------------------------------------------!
subroutine bcs_l_block_matvec(A,x,y)                                       !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)

    call A%l_block_matvec_impl(x,y)

end subroutine bcs_l_block_matvec



!--------------------------------------------------------------------------!
subroutine bcs_l_block_matvec_t(A,x,y)                                     !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)

    call A%l_block_matvec_t_impl(x,y)

end subroutine bcs_l_block_matvec_t



!--------------------------------------------------------------------------!
subroutine bcs_r_block_matvec(A,x,y)                                       !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)

    call A%r_block_matvec_impl(x,y)

end subroutine bcs_r_block_matvec



!--------------------------------------------------------------------------!
subroutine bcs_r_block_matvec_t(A,x,y)                                     !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)

    call A%r_block_matvec_t_impl(x,y)

end subroutine bcs_r_block_matvec_t



!--------------------------------------------------------------------------!
subroutine bcs_set_block_not_preallocated(A,i,j)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: i,j
!    real(dp), intent(in) :: val
    ! local variables
    real(dp), allocatable :: val_temp(:,:,:)
    integer :: k,n1,n2

    if (A%orientation=='col') then
        call A%g%add_edge(j,i)
    else
        call A%g%add_edge(i,j)
    endif

    k = A%find_block(i,j)

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    allocate(val_temp(n1,n2,A%nnz+1))

    val_temp(:,:,1:k-1) = A%val(:,:,1:k-1)
    val_temp(:,:,k) = 0.0_dp
    val_temp(:,:,k+1:A%nnz+1) = A%val(:,:,k:A%nnz)
    call move_alloc(from=val_temp, to=A%val)

    A%nnz = A%nnz+1
    A%max_degree = A%g%max_degree

end subroutine bcs_set_block_not_preallocated





!==========================================================================!
!==========================================================================!
!==== Implementations of matrix operations                             ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
function bcsr_find_block(A,i,j)                                            !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: bcsr_find_block

    bcsr_find_block = A%g%find_edge(i,j)

end function bcsr_find_block



!--------------------------------------------------------------------------!
function bcsc_find_block(A,i,j)                                            !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(in) :: A
    integer, intent(in) :: i,j
    integer :: bcsc_find_block

    bcsc_find_block = A%g%find_edge(j,i)

end function bcsc_find_block



!--------------------------------------------------------------------------!
subroutine bcs_matrix_permute_ptrs(A,p)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i,k,ptr(A%g%n+1)
    real(dp) :: val(size(A%val,1),size(A%val,2),A%nnz)

    do i=1,A%g%n
        ptr(p(i)+1) = A%g%ptr(i+1)-A%g%ptr(i)
    enddo

    ptr(1) = 1
    do i=1,A%g%n
        ptr(i+1) = ptr(i+1)+ptr(i)
    enddo

    do i=1,A%g%n
        do k=0,A%g%ptr(i+1)-A%g%ptr(i)-1
            val(:,:,ptr(p(i)+k)) = A%val(:,:,A%g%ptr(i)+k)
        enddo
    enddo

    A%val = val

    call A%g%left_permute(p)

end subroutine bcs_matrix_permute_ptrs



!--------------------------------------------------------------------------!
subroutine bcs_matrix_permute_vals(A,p)                                    !
!--------------------------------------------------------------------------!
    class(bcs_matrix), intent(inout) :: A
    integer, intent(in) :: p(:)

    call A%g%right_permute(p)

end subroutine bcs_matrix_permute_vals



!--------------------------------------------------------------------------!
subroutine bcsr_matvec(A,x,y)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(size(A%val,1))

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
!            do l=1,A%nc
!                z = z+A%val(:,l,k)*x(A%nc*(j-1)+l)
!            enddo
            z = z+matmul(A%val(:,:,k),x(n2*(j-1)+1:n2*j))
        enddo
        y(n1*(i-1)+1:n1*i) = z
    enddo

end subroutine bcsr_matvec



!--------------------------------------------------------------------------!
subroutine bcsc_matvec(A,x,y)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(size(A%val,1))

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    do j=1,A%g%n
        z = x(n1*(j-1)+1:n1*j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
            y(n2*(i-1)+1:n2*i) = y(n2*(i-1)+1:n2*i)+matmul(z,A%val(:,:,k))
!            do l=1,n1
!                y(A%nc*(j-1)+l) = y(A%nc*(j-1)+1) &
!                    & +dot_product(A%val(:,k,l),z)
!            enddo
        enddo
    enddo

end subroutine bcsc_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_block_matvec(A,x,y)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(size(A%val,1))

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
            z = z+matmul(A%val(:,:,k),x(:,j))
        enddo
        y(:,i) = z
    enddo

end subroutine bcsr_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsc_block_matvec(A,x,y)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l
    real(dp) :: z(size(A%val,1))

    y = 0.0_dp
    do j=1,A%g%n
        z = x(:,j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
            y(:,i) = y(:,i)+matmul(z,A%val(:,:,k))
        enddo
    enddo

end subroutine bcsc_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_l_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(A%nr)

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
!            do l=1,A%nc
!                z = z+A%val(:,l,k)*x(A%nc*(j-1)+l)
!            enddo
            z = z+matmul(A%val(:,:,k),x(n2*(j-1)+1:n2*j))
        enddo
        y(:,i) = z
    enddo

end subroutine bcsr_l_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsc_l_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:,:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(size(A%val,1))

    do j=1,A%g%n
        z = x(n1*(j-1)+1:n1*j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
!            do l=1,A%nc
!                y(l,j) = y(l,j)+dot_product(A%val(:,l,k),z)
!            enddo
            y(:,i) = y(:,i)+matmul(z,A%val(:,:,k))
        enddo
    enddo

end subroutine bcsc_l_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsr_r_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(A%nr)

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    do i=1,A%g%n
        z = 0.0_dp
        do k=A%g%ptr(i),A%g%ptr(i+1)-1
            j = A%g%node(k)
!            do l=1,A%nc
!                z = z+A%val(:,l,k)*x(l,j)
!            enddo
            z = z+matmul(A%val(:,:,k),x(:,j))
        enddo
        y(n1*(i-1)+1:n1*i) = z
    enddo

end subroutine bcsr_r_block_matvec



!--------------------------------------------------------------------------!
subroutine bcsc_r_block_matvec(A,x,y)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(bcs_matrix), intent(in) :: A
    real(dp), intent(in)  :: x(:,:)
    real(dp), intent(out) :: y(:)
    ! local variables
    integer :: i,j,k,l,n1,n2
    real(dp) :: z(size(A%val,1))

    n1 = size(A%val,1)
    n2 = size(A%val,2)

    do j=1,A%g%n
        z = x(:,j)
        do k=A%g%ptr(j),A%g%ptr(j+1)-1
            i = A%g%node(k)
!            do l=1,A%nc
!                y(A%nc*(j-1)+l) = y(A%nc*(j-1)+l) &
!                    & +dot_product(A%val(:,l,k),z)
!            enddo
            y(n2*(i-1)+1:n2*i) = y(n2*(i-1)+1:n2*i)+matmul(z,A%val(:,:,k))
        enddo
    enddo

end subroutine bcsc_r_block_matvec






end module bcs_matrices
