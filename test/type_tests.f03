program type_tests

use types

implicit none

    type(linked_list) :: list
    integer :: i,j

    print *, 'Setting list = {1, ..., 10}'
    do i=1,10
        call list%append(i)
    enddo

    print *, 'Getting entry #5'
    i = list%get_entry(5)
    if (i/=5) then
        print *, 'list(5) should be = 5;'
        print *, 'value found: ',i
    endif

    print *, 'Deleting entry #5'
    call list%delete_entry(5)
    i = list%get_value(5)
    if (i/=6) then
        print *, 'list(5) should now be = 5;'
        print *, 'values found: ',i
    endif


end program type_tests
