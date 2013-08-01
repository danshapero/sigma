program utility_test

use utility_algorithms

implicit none



    integer :: list(10),order(10)

    list = [10, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    print *, order(list)
    print *, list(order(list))




end program utility_test
