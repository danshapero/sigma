program util_tests

use util

implicit none



    integer :: i, list(10), p(10)
    real(dp) :: det, A(5,5)
    

    ! List sorting
    list = [10, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    p = order(list)

    do i = 1, 9
        if (list(p(i)) > list(p(i+1))) then
            print *, "Sorting list failed!"
            call exit(1)
        endif
    enddo


    ! Finding elements in sorted list
    list = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100]
    if (member(list, 3)) then
        print *, "Finding element in list returned false positive!"
        call exit(1)
    endif

    do i = 1, 10
        if (.not. member(list, i**2)) then
            print *, "Finding element in list failed to find element", i**2
            call exit(1)
        endif
    enddo


    ! Determinant
    A = 0.0_dp
    do i = 1, 4
        A(i, i+1) = -1.0_dp
        A(i+1, i) = -1.0_dp
        A(i, i) = 2.0_dp
    enddo
    A(5,5) = 2.0_dp

    det = determinant(A)
    if (dabs(det - 6.0_dp) > 1.0e-16) then
        print *, 'Determinant of A should be = 6.0'
        print *, 'Value found: ', det
        call exit(1)
    endif


end program util_tests
