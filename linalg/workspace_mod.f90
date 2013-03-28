!--------------------------------------------------------------------------!
module workspace_mod                                                       !
!--------------------------------------------------------------------------!

    implicit none


type :: workspace
    integer :: nn,nv
    real(kind(1d0)), pointer :: w(:,:)
    logical, allocatable :: in_use(:)
contains
    procedure :: init_workspace
    procedure :: find_free_work_vectors
    procedure :: point_to_work_vector
    procedure :: free_work_vectors
end type workspace



contains



!--------------------------------------------------------------------------!
subroutine init_workspace(work,nn,nv)                                      !
!--------------------------------------------------------------------------!
    implicit none
    class(workspace), intent(inout) :: work
    integer, intent(in) :: nn,nv

    work%nn = nn
    work%nv = nv

    allocate( work%w(nn,nv), work%in_use(nv) )
    work%in_use = .false.

end subroutine init_workspace



!--------------------------------------------------------------------------!
integer function find_free_work_vectors(work,k)                            !
!--------------------------------------------------------------------------!
    implicit none
    class(workspace), intent(in) :: work
    integer, intent(in) :: k
    ! local variables
    integer :: i,j
    logical :: free

    do i=1,work%nv-k+1
        free = .true.
        do j=0,k-1
            free = free .and. .not.work%in_use(i+j)
        enddo

        if (free) then
            find_free_work_vectors = i
            return
        endif
    enddo

    find_free_work_vectors = -1
    return
    
end function find_free_work_vectors



!--------------------------------------------------------------------------!
subroutine point_to_work_vector(work,ind,ptr)                              !
!--------------------------------------------------------------------------!
    implicit none
    class(workspace), intent(inout) :: work
    integer, intent(in) :: ind
    real(kind(1d0)), pointer, intent(inout) :: ptr(:)

    ptr => work%w(:,ind)
    work%in_use(ind) = .true.

end subroutine point_to_work_vector



!--------------------------------------------------------------------------!
subroutine free_work_vectors(work,ind,k)                                   !
!--------------------------------------------------------------------------!
    implicit none
    class(workspace), intent(inout) :: work
    integer, intent(in) :: ind,k
    ! local variables
    integer :: i

    do i=ind,ind+k-1
        work%in_use(i) = .false.
    enddo

end subroutine free_work_vectors



end module workspace_mod
