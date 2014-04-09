module types

implicit none

integer, parameter :: dp=kind(0.d0)
real(dp), parameter :: pi=3.1415926535897932_dp



!--------------------------------------------------------------------------!
type :: dynamic_array                                                      !
!--------------------------------------------------------------------------!
    integer :: length, capacity, min_capacity
    integer, allocatable :: array(:)
contains
    procedure :: init => dynamic_array_init
    procedure :: get_entry => dynamic_array_get_entry
    procedure :: set_entry => dynamic_array_set_entry
    procedure :: push => dynamic_array_push
    procedure :: pop => dynamic_array_pop
    procedure :: peek => dynamic_array_peek
    procedure :: free => dynamic_array_free
end type dynamic_array



!--------------------------------------------------------------------------!
type, extends(dynamic_array) :: circular_array                             !
!--------------------------------------------------------------------------!
    integer :: start
contains
    procedure :: init => circular_array_init
    procedure :: get_entry => circular_array_get_entry
    procedure :: set_entry => circular_array_set_entry
    ! stack-like operations
    procedure :: push => circular_array_push
    procedure :: pop => circular_array_pop
    procedure :: peek => circular_array_peek
    ! queue-like operations
    procedure :: enqueue => circular_array_enqueue
    procedure :: dequeue => circular_array_dequeue
    procedure :: front => circular_array_front
    ! auxiliary operations
    procedure, private :: circular_array_expand
    procedure, private :: circular_array_contract
end type circular_array



contains




!==========================================================================!
!==========================================================================!
!==== Methods for dynamic arrays                                       ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
subroutine dynamic_array_init(a,capacity,min_capacity)                     !
!--------------------------------------------------------------------------!
    class(dynamic_array), intent(inout) :: a
    integer, intent(in), optional :: capacity, min_capacity

    a%length = 0
    a%min_capacity = 4
    a%capacity = 32

    if (present(capacity)) then
        a%capacity = capacity
    endif

    if (present(min_capacity)) then
        a%min_capacity = min_capacity
    endif

    allocate(a%array(a%capacity))
    a%array = 0

end subroutine dynamic_array_init



!--------------------------------------------------------------------------!
elemental function dynamic_array_get_entry(a,i) result(val)                !
!--------------------------------------------------------------------------!
    class(dynamic_array), intent(in) :: a
    integer, intent(in) :: i
    integer :: val

    val = a%array(i)

end function dynamic_array_get_entry



!--------------------------------------------------------------------------!
elemental subroutine dynamic_array_set_entry(a,i,val)                      !
!--------------------------------------------------------------------------!
    class(dynamic_array), intent(inout) :: a
    integer, intent(in) :: i,val

    a%array(i) = val

end subroutine dynamic_array_set_entry



!--------------------------------------------------------------------------!
subroutine dynamic_array_push(a,val)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(dynamic_array), intent(inout) :: a
    integer, intent(in) :: val
    ! local variables
    integer, allocatable :: array(:)

    if (a%length==a%capacity) then
        allocate(array(2*a%capacity))
        a%capacity = 2*a%capacity
        array = 0
        array(1:a%length) = a%array(1:a%length)
        call move_alloc(from=array, to=a%array)
    endif

    a%length = a%length+1
    a%array(a%length) = val

end subroutine dynamic_array_push



!--------------------------------------------------------------------------!
function dynamic_array_pop(a) result(val)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(dynamic_array), intent(inout) :: a
    integer :: val
    ! local variables
    integer, allocatable :: array(:)

    if (a%length<a%capacity/4 .and. a%capacity/2>=a%min_capacity) then
        allocate(array(a%capacity/2))
        a%capacity = a%capacity/2
        array(1:a%length) = a%array(1:a%length)
        call move_alloc(from=array, to=a%array)
    endif

    if (a%length>0) then
        val = a%array(a%length)

        a%array(a%length) = 0
        a%length = a%length-1
    else
        print *, 'Error: cannot pop from an empty dynamic array'
        call exit(1)
    endif

end function dynamic_array_pop



!--------------------------------------------------------------------------!
function dynamic_array_peek(a) result(val)                                 !
!--------------------------------------------------------------------------!
    class(dynamic_array), intent(in) :: a
    integer :: val

    val = 0
    if (a%length>0) val = a%array(a%length)

end function dynamic_array_peek



!--------------------------------------------------------------------------!
subroutine dynamic_array_free(a)                                           !
!--------------------------------------------------------------------------!
    class(dynamic_array), intent(inout) :: a

    deallocate(a%array)

    a%length = 0
    a%capacity = 0
    a%min_capacity = 0

end subroutine dynamic_array_free




!==========================================================================!
!==========================================================================!
!==== Methods for circular arrays                                      ====!
!==========================================================================!
!==========================================================================!



!--------------------------------------------------------------------------!
subroutine circular_array_init(a,capacity,min_capacity)                    !
!--------------------------------------------------------------------------!
     class(circular_array), intent(inout) :: a
     integer, intent(in), optional :: capacity, min_capacity

     call a%dynamic_array%init(capacity,min_capacity)

     a%start = a%capacity/2+1

end subroutine circular_array_init



!--------------------------------------------------------------------------!
elemental function circular_array_get_entry(a,i) result(val)               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(in) :: a
    integer, intent(in) :: i
    integer :: val
    ! local variables
    integer :: k

    ! Find the offset of entry i from the start of the array
    k = a%start+i-1

    ! Find the actual index in a%array where this entry would be stored
    k = mod(k-1,a%capacity)+1

    ! Return the value stored at that index
    val = a%array(k)

end function circular_array_get_entry



!--------------------------------------------------------------------------!
elemental subroutine circular_array_set_entry(a,i,val)                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    integer, intent(in) :: i, val
    ! local variables
    integer :: k

    k = a%start+i-1
    k = mod(k-1,a%capacity)+1
    a%array(k) = val

end subroutine circular_array_set_entry



!--------------------------------------------------------------------------!
subroutine circular_array_push(a,val)                                      !
!--------------------------------------------------------------------------!
     ! input/output variables
     class(circular_array), intent(inout) :: a
     integer, intent(in) :: val
     ! local variables
     integer :: k

     associate(start => a%start, length => a%length)

     ! If the array is at capacity, expand the array
     if (length==a%capacity) call a%circular_array_expand()

     ! Increment the length of the array
     length = length+1

     ! Find the index in which to put the new element
     k = start+length-1
     k = mod(k-1,a%capacity)+1

     ! Store the new value in the array
     a%array(k) = val

     end associate

end subroutine circular_array_push



!--------------------------------------------------------------------------!
function circular_array_pop(a) result(val)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    integer :: val
    ! local variables
    integer :: k

    associate(start => a%start, length => a%length)

    ! If the array is using less than 1/4 of its capacity, contract it to
    ! half its size
    if (length<a%capacity/4 .and. a%capacity/2>=a%min_capacity) then
        call a%circular_array_contract()
    endif

    if (length>0) then
        k = start+length-1
        k = mod(k-1,a%capacity)+1

        val = a%array(k)

        a%array(k) = 0
        length = length-1
    else
        print *, 'Error: cannot pop from an empty circular array'
        call exit(1)
    endif

    end associate

end function circular_array_pop



!--------------------------------------------------------------------------!
function circular_array_peek(a) result(val)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(in) :: a
    integer :: val
    ! local variables
    integer :: k

    k = a%start+a%length-1
    k = mod(k-1,a%capacity)+1

    ! print *, 'Hi!'

    val = 0
    if (a%length>0) val = a%array(k)

end function circular_array_peek



!--------------------------------------------------------------------------!
subroutine circular_array_enqueue(a,val)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    integer, intent(in) :: val
    ! local variables
    integer :: k

    associate(start => a%start, length => a%length)

    ! If the array is at capacity, expand the array
    if (length==a%capacity) call a%circular_array_expand()

    ! Find the index in which to put the new element
    k = start+a%capacity-1
    k = mod(k-1,a%capacity)+1

    ! Store the new value in the array
    a%array(k) = val

    ! Increment the length of the array and move the new start position
    ! either down by one or wrap it around the end of the underlying linear
    ! array
    start = k
    length = length+1

    end associate

end subroutine circular_array_enqueue



!--------------------------------------------------------------------------!
function circular_array_dequeue(a) result(val)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    integer :: val
    ! local variables
    integer :: k

    associate(start => a%start, length => a%length)

    ! If the array is using less than 1/4 of its capacity, contract it to
    ! half its size
    if (length<a%capacity/4 .and. a%capacity/2>=a%min_capacity) then
        call a%circular_array_contract()
    endif

    if (length>0) then
        ! The value to be returned is the front of the circular array
        val = a%array(start)

        ! Clear the front of the array and decrement the length
        a%array(start) = 0
        length = length-1

        ! Increment the start, modulo the array's capacity
        k = start+1
        k = mod(k-1,a%capacity)+1

        start = k
    else
        print *, 'Error: cannot dequeue from an empty circular array'
        call exit(1)
    endif

    end associate

end function circular_array_dequeue



!--------------------------------------------------------------------------!
function circular_array_front(a) result(val)                               !
!--------------------------------------------------------------------------!
    class(circular_array), intent(in) :: a
    integer :: val

    val = 0
    if (a%length>0) val = a%array(a%start)

end function circular_array_front



!--------------------------------------------------------------------------!
subroutine circular_array_expand(a)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    ! local variables
    integer :: last,num_before_end,num_after_end
    integer, allocatable :: array(:)

    associate(start => a%start, length => a%length)

    allocate(array(2*a%capacity))
    array = 0

    last = start+length-1
    num_before_end = min(length,a%capacity-start+1)
    num_after_end  = length-num_before_end

    array(1:num_before_end) = a%array(start:start+num_before_end-1)
    array(num_before_end+1:length) = a%array(1:num_after_end)
    start = 1

    call move_alloc(from=array, to=a%array)
    a%capacity = 2*a%capacity

    end associate

end subroutine circular_array_expand



!--------------------------------------------------------------------------!
subroutine circular_array_contract(a)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(circular_array), intent(inout) :: a
    ! local variables
    integer :: last,num_before_end,num_after_end
    integer, allocatable :: array(:)

    associate(start => a%start, length => a%length)

    allocate(array(a%capacity/2))
    array = 0

    last = start+length-1
    num_before_end = min(length,a%capacity-start+1)
    num_after_end  = length-num_before_end

    array(1:num_before_end) = a%array(start:start+num_before_end-1)
    array(num_before_end+1:length) = a%array(1:num_after_end)
    start = 1

    call move_alloc(from=array, to=a%array)
    a%capacity = a%capacity/2

    end associate

end subroutine circular_array_contract




end module types
