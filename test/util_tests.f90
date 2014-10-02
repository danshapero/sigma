program util_tests

use util

implicit none



    integer :: i,list(10),p(10)
    real(dp) :: det, A(5,5)
    

    list = [10, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    p = order(list)
    !print *, list(order(list))


    A = 0.0_dp
    do i=1,4
        A(i,i+1) = -1.0_dp
        A(i+1,i) = -1.0_dp
        A(i,i) = 2.0_dp
    enddo
    A(5,5) = 2.0_dp
    det = determinant(A)
    if (dabs(det - 6.0_dp) > 1.0e-16) then
        print *, 'Determinant of A should be = 6.0'
        print *, 'Value found: ',det
    endif


end program util_tests
