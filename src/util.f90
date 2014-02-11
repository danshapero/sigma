module util

use types

implicit none

contains


!==========================================================================!
!==========================================================================!
!==== Assorted useful algorithms                                       ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
function order(list)                                                       !
!--------------------------------------------------------------------------!
! Produce a permutation of a list which puts it in order:                  !
!       list( order(list) ) = the sorted array list                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    integer, intent(in) :: list(:)
    integer :: order(size(list))
    ! local variables
    integer :: i,j,key,indx

    do i=1,size(list)
        order(i) = i
    enddo

    do i=2,size(list)
        indx = order(i)
        key = list(indx)
        do j=i-1,1,-1
            if (list(order(j))<=key) exit
            order(j+1) = order(j)
        enddo
        order(j+1) = indx
    enddo

end function order




!--------------------------------------------------------------------------!
function determinant(A)                                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    real(dp), intent(in) :: A(:,:)
    real(dp) :: determinant
    ! local variables
    integer :: i, ierr, ipiv(size(A,1))
    real(dp) :: B(size(A,1),size(A,2))

    B = A
    call dgetrf(size(A,1),size(A,2),B,size(A,1),ipiv,ierr)

    determinant = 1.0_dp
    do i=1,size(A,2)
        determinant = determinant*B(i,i)
        determinant = determinant*(-1)**(i+ipiv(i))
    enddo

end function determinant



!--------------------------------------------------------------------------!
subroutine init_seed()                                                     !
!--------------------------------------------------------------------------!
! Seed a pseudo-random number generator using the current date and time    !
! Algorithm courtesy of the Fortran wiki, see                              !
!     http://fortranwiki.org/fortran/show/init_seed                        !
!--------------------------------------------------------------------------!
    integer :: n, ival(8), v(3), i
    integer, allocatable :: seed(:)

    call date_and_time(values=ival)

    

    v(1) = ival(8) + 2048*ival(7)
    v(2) = ival(6) + 64*ival(5)     ! value(4) isn't really 'random'
    v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)

    call random_seed(size=n)
    allocate(seed(n))
    call random_seed()   ! Give the seed an implementation-dependent kick
    call random_seed(get=seed)

    do i=1,n
        seed(i) = seed(i) + v(mod(i-1, 3) + 1)
    enddo

    call random_seed(put=seed)

    deallocate(seed)

end subroutine








end module util
