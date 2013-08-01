module util

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
    implicit none
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









end module util
