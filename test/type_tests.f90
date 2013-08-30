program type_tests

use types

implicit none

    type(linked_list) :: list
    integer :: i,j

    do i=1,10
        call list%append(i)
    enddo

    i = list%get_entry(5)
    if (i/=5) then
        print *, 'list(5) should be = 5;'
        print *, 'value found: ',i
    endif

    call list%delete_entry(5)
    i = list%get_value(5)
    if (i/=6) then
        print *, 'list(5) should now be = 5;'
        print *, 'values found: ',i
    endif


end program type_tests
