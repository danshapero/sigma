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








end module util
