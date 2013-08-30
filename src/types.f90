module types

implicit none

integer, parameter :: dp=kind(0.d0)
real(dp), parameter :: pi=3.1415926535897932_dp



type :: node
    integer :: val
    type(node), pointer :: next => null()
end type node

type :: linked_list
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


contains



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



end module types
