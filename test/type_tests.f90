program type_tests

use types

implicit none

    type(dynamic_array) :: a
    type(circular_array) :: c
    integer :: i,j


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

    call a%free()


    !----------------------------------------------------------------------!
    ! Circular array tests                                                 !
    !----------------------------------------------------------------------!
    call c%init()
    if (c%length/=0) then
        print *, 'Length of just-initialized circular array should be 0'
        call exit(1)
    endif

    call c%push(42)
    if (c%length/=1) then
        print *, 'Length of circular array after pushing a single int '
        print *, 'should be 1, we found: ', c%length
        call exit(1)
    endif
    call c%enqueue(88)

    i = c%front()
    j = c%peek()
    if (i/=88 .or. j/=42) then
        print *, 'We pushed 42 and enqueued 88; looking at the front and '
        print *, 'the end of the circular array should yield 42 and 88, '
        print *, 'instead we got: '
        print *, 'Front:',i
        print *, 'Back:',j
        call exit(1)
    endif

    do i=1,18
        call c%push(i)
        call c%enqueue(i)
    enddo

    do i=18,1,-1
        j = c%dequeue()
        if (j/=i) then
            print *, 'Should have dequeued ',i,'instead got',j
            call exit(1)
        endif
    enddo

    do i=18,1,-1
        j = c%pop()
    enddo

    call c%free()


    !----------------------------------------------------------------------!
    ! End of tests                                                         !
    !----------------------------------------------------------------------!
    call exit(0)


end program type_tests
