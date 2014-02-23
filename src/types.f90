module types

implicit none

integer, parameter :: dp=kind(0.d0)
real(dp), parameter :: pi=3.1415926535897932_dp


type :: node
    integer :: val
    type(node), pointer :: next => null()
end type node


!--------------------------------------------------------------------------!
type :: linked_list                                                        !
!--------------------------------------------------------------------------!
    integer :: length = 0
    type(node), pointer :: head => null()
contains
    procedure :: append => linked_list_append
    procedure :: prepend => linked_list_prepend
    procedure :: set_value => linked_list_set_value
    procedure :: delete_entry => linked_list_delete_entry
    procedure :: delete_value => linked_list_delete_value
    procedure :: get_entry => linked_list_get_entry
    procedure :: get_value => linked_list_get_value
    procedure :: free => linked_list_free
end type linked_list


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
!==== Methods for linked lists                                         ====!
!==========================================================================!
!==========================================================================!


!--------------------------------------------------------------------------!
function new_node(i)                                                       !
!--------------------------------------------------------------------------!
    integer, intent(in) :: i
    type(node), pointer :: new_node

    allocate(new_node)
    new_node%val = i
    nullify(new_node%next)

end function new_node



!--------------------------------------------------------------------------!
subroutine linked_list_append(list,i)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(inout) :: list
    integer, intent(in) :: i
    ! local variables
    type(node), pointer :: current

    if (.not.associated(list%head)) then
        list%head => new_node(i)
        list%length = 1
    else
        current => list%head
        do while(associated(current%next))
            current => current%next
        enddo
        current%next => new_node(i)
        list%length = list%length+1
    endif

end subroutine linked_list_append



!--------------------------------------------------------------------------!
subroutine linked_list_prepend(list,i)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(inout) :: list
    integer, intent(in) :: i
    ! local variables
    type(node), pointer :: current

    if (.not.associated(list%head)) then
        list%head => new_node(i)
        list%length = 1
    else
        current => new_node(i)
        current%next => list%head
        list%head => current
        list%length = list%length+1
    endif

end subroutine linked_list_prepend



!--------------------------------------------------------------------------!
subroutine linked_list_set_value(list,i,val)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(inout) :: list
    integer, intent(in) :: i, val
    ! local variables
    integer :: k
    type(node), pointer :: current

    if (.not.associated(list%head)) then
        print *, 'Linked list is empty'
        call exit(1)
    endif
    current => list%head
    do k=1,i-1
        if (.not.associated(current%next)) then
            print *, 'Linked list index ',i,' out of bounds'
            call exit(1)
        endif
        current => current%next
    enddo
    current%val = val

end subroutine linked_list_set_value



!--------------------------------------------------------------------------!
subroutine linked_list_delete_entry(list,i)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(inout) :: list
    integer, intent(in) :: i
    ! local variables
    integer :: k
    type(node), pointer :: previous,current

    ! Check if the list has any elements
    if (associated(list%head)) then
        ! If we're deleting the first element,
        if (i==1) then
            ! check that there's a second element to make into the new head.
            if (associated(list%head%next)) then
                current => list%head%next
                deallocate(list%head)
                list%head => current
            ! Otherwise, delete the head of the list.
            else
                deallocate(list%head)
            endif
        ! If we're not deleting the first element,
        else
            previous => list%head
            current => list%head
            do k=1,i-1
                previous => current
                current => current%next
            enddo
            previous%next => current%next
            deallocate(current)
        endif

        list%length = list%length-1
    endif

end subroutine linked_list_delete_entry



!--------------------------------------------------------------------------!
subroutine linked_list_delete_value(list,i)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(inout) :: list
    integer, intent(in) :: i
    ! local variables
    type(node), pointer :: previous, current

    if (associated(list%head)) then
        if (list%head%val==i) then
            if (associated(list%head%next)) then
                current => list%head%next
                deallocate(list%head)
                list%head => current
            else
                deallocate(list%head)
            endif
            list%length = list%length-1
        else
            previous => list%head
            current => list%head
            do while(associated(current%next))
                previous => current
                current => current%next
                if (current%val==i) then
                    previous%next => current%next
                    deallocate(current)
                    list%length = list%length-1
                    exit
                endif
            enddo
        endif
    endif

end subroutine linked_list_delete_value



!--------------------------------------------------------------------------!
function linked_list_get_entry(list,val)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(in) :: list
    integer, intent(in) :: val
    integer :: linked_list_get_entry
    ! local variables
    integer :: k
    type(node), pointer :: current

    linked_list_get_entry = -1
    current => list%head
    k = 1
    do while(associated(current%next))
        if (current%val==val) then
            linked_list_get_entry = k
            exit
        endif
        current => current%next
        k = k+1
    enddo

end function linked_list_get_entry


!--------------------------------------------------------------------------!
function linked_list_get_value(list,i)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(linked_list), intent(in) :: list
    integer, intent(in) :: i
    integer :: linked_list_get_value
    ! local variables
    integer :: k
    type(node), pointer :: current

    if (.not.associated(list%head)) then
        print *, 'Linked list is empty'
        call exit(1)
    endif
    current => list%head
    do k=1,i-1
        if (.not.associated(current%next)) then
            print *, 'Linked list index ',i,'out of bounds'
            call exit(1)
        endif
        current => current%next
    enddo

    linked_list_get_value = current%val

end function linked_list_get_value



!--------------------------------------------------------------------------!
subroutine linked_list_free(list)                                          !
!--------------------------------------------------------------------------!
    class(linked_list), intent(inout) :: list
    integer :: k,length

    length = list%length

    do k=1,length
        call list%delete_entry(1)
    enddo

end subroutine linked_list_free




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
        a%capacity = max(capacity,a%capacity)
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
