program type_tests

use types

implicit none

    type(linked_list) :: list
    type(dynamic_array) :: a
    integer :: i,j


    !----------------------------------------------------------------------!
    ! Linked list tests                                                    !
    !----------------------------------------------------------------------!
    do i=1,10
        call list%append(i)
    enddo

    i = list%get_entry(5)
    if (i/=5) then
        print *, 'list(5) should be = 5;'
        print *, 'value found: ',i
        call exit(1)
    endif

    call list%delete_entry(5)
    i = list%get_value(5)
    if (i/=6) then
        print *, 'list(5) should now be = 5;'
        print *, 'values found: ',i
        call exit(1)
    endif


    !----------------------------------------------------------------------!
    ! Dynamic array tests                                                  !
    !----------------------------------------------------------------------!
    call a%init()
    if (a%length/=0) then
        print *, 'Length of just-initialized dynamic array should be 0'
        call exit(1)
    endif

    call a%push(1)
    if (a%length/=1) then
        print *, 'Length of dynamic array after pushing a single int '
        print *, 'should be 1, we found: ', a%length
        call exit(1)
    endif

    i = a%pop()
    if (i/=1) then
        print *, 'Should have popped 1 off the dynamic array; found',i
        call exit(1)
    endif

    do i=1,32
        call a%push(i)
        j = a%peek()
        if (j/=i) then
            print *, 'Just pushed',i
            print *, 'Peeking at the dynamic array should yield',i
            print *, 'Instead got',j
            call exit(1)
        endif

        if (a%length/=i) then
            print *, 'Length not properly updated'
            call exit(1)
        endif
    enddo


end program type_tests
