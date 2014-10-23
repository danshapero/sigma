module coo_graphs

use graph_interfaces
use types, only: dynamic_array

implicit none



!--------------------------------------------------------------------------!
type, extends(graph_interface) :: coo_graph                                !
!--------------------------------------------------------------------------!
    type(dynamic_array) :: edges(2)
    integer, private :: ne
contains
    !--------------
    ! Constructors
    procedure :: init_empty => coo_graph_init
    procedure :: init_from_iterator => coo_graph_init_from_iterator

    !-----------
    ! Accessors
    procedure :: get_num_edges => coo_get_num_edges
    procedure :: get_degree => coo_get_degree
    procedure :: get_max_degree => coo_get_max_degree
    procedure :: get_neighbors => coo_get_neighbors
    procedure :: connected => coo_connected
    procedure :: find_edge => coo_find_edge

    !---------------
    ! Edge iterator
    procedure :: make_cursor => coo_make_cursor
    procedure :: get_edges => coo_get_edges

    !----------
    ! Mutators
    procedure :: add_edge => coo_add_edge
    procedure :: delete_edge => coo_delete_edge
    procedure :: left_permute => coo_graph_left_permute
    procedure :: right_permute => coo_graph_right_permute

    !-------------
    ! Destructors
    procedure :: destroy => coo_destroy

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => coo_dump_edges
end type coo_graph





contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine coo_graph_init(g, n, m)                                         !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m

    call g%add_reference()

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    call g%edges(1)%init(capacity = 4)
    call g%edges(2)%init(capacity = 4)

    ! At initialization, the number of edges and max degree is zero
    g%ne = 0

end subroutine coo_graph_init



!--------------------------------------------------------------------------!
subroutine coo_graph_init_from_iterator(g, n, m, &                         !
                                            & iterator, get_cursor, trans) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: n, m
    procedure(edge_iterator) :: iterator
    procedure(new_cursor) :: get_cursor
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ord(2)
    integer :: i, j, k
    integer :: num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    call g%add_reference()

    ord = [1, 2]
    if (present(trans)) then
        if (trans) ord = [2, 1]
    endif

    g%n = n
    g%m = m

    cursor = get_cursor()

    g%ne = cursor%last

    call g%edges(1)%init(capacity = g%ne + 16)
    call g%edges(2)%init(capacity = g%ne + 16)

    do while (.not. cursor%done())
        call iterator(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(ord(1), k)
            j = edges(ord(2), k)

            call g%edges(1)%push(i)
            call g%edges(2)%push(j)
        enddo
    enddo


end subroutine coo_graph_init_from_iterator




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function coo_get_num_edges(g) result(k)                                    !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    integer :: k

    k = g%ne

end function coo_get_num_edges



!--------------------------------------------------------------------------!
function coo_get_degree(g, i) result(d)                                    !
!--------------------------------------------------------------------------!
!    NOTE: This is an extremely inefficient procedure. In order to be able !
! have a COO graph require only O(ne) storage, and be able to add new      !
! edges in O(1) time, we have to allow for duplicate edges. That makes     !
! checking the degree of a vertex require O(ne) time.                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d
    ! local variables
    integer :: j, k, l
    logical :: found
    type(dynamic_array) :: neighbors

    call neighbors%init(capacity = 4, min_capacity = 2)

    ! Loop through all the edges of `g`
    do k = 1, g%ne
        ! If the current edge starts at `i`,
        if (g%edges(1)%get_entry(k) == i) then
            ! store the ending vertex `j` of that edge.
            j = g%edges(2)%get_entry(k)

            ! Check to make sure that vertex `j` has not already been
            ! entered into the list of neighbors.
            found = .false.

            do l = 1, neighbors%length
                if (neighbors%get_entry(l) == j) found = .true.
            enddo

            if (.not. found) call neighbors%push(j)
        endif
    enddo

    ! Return the length of the neighbors array
    d = neighbors%length

end function coo_get_degree



!--------------------------------------------------------------------------!
function coo_get_max_degree(g) result(k)                                    !
!--------------------------------------------------------------------------!
!     NOTE: This procedure is very inefficient, being an n-fold repetition !
! of the already very inefficient `degree` method.                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer :: k
    ! local variables
    integer :: i

    k = 0

    do i = 1, g%n
        k = max(k, g%get_degree(i))
    enddo

end function coo_get_max_degree



!--------------------------------------------------------------------------!
subroutine coo_get_neighbors(g, neighbors, i)                              !
!--------------------------------------------------------------------------!
!    See comment for `coo_degree` procedure; this procedure is very, very  !
! inefficient.                                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: j, k, l
    logical :: found
    type(dynamic_array) :: nbrs

    call nbrs%init(capacity = 4, min_capacity = 2)

    do k = 1, g%ne
        if (g%edges(1)%get_entry(k) == i) then
            j = g%edges(2)%get_entry(k)

            found = .false.
            do l = 1, nbrs%length
                if (nbrs%get_entry(l) == j) found = .true.
            enddo

            if (.not. found) call nbrs%push(j)
        endif
    enddo

    neighbors = 0

    do k = 1, nbrs%length
        neighbors(k) = nbrs%get_entry(k)
    enddo

end subroutine coo_get_neighbors



!--------------------------------------------------------------------------!
function coo_connected(g, i, j)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i, j
    logical :: coo_connected
    ! local variables
    integer :: k

    coo_connected = .false.

    do k = 1,g%ne
        if (g%edges(1)%get_entry(k) == i &        
                                & .and. g%edges(2)%get_entry(k) == j) then
            coo_connected = .true.
        endif
    enddo

end function coo_connected



!--------------------------------------------------------------------------!
function coo_find_edge(g, i, j)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: i, j
    integer :: coo_find_edge
    ! local variables
    integer :: k

    coo_find_edge = -1

    do k = 1, g%ne
        if (g%edges(1)%get_entry(k) == i &
                                & .and. g%edges(2)%get_entry(k) == j) then
            coo_find_edge = k
        endif
    enddo

end function coo_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function coo_make_cursor(g) result(cursor)                                 !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    type(graph_edge_cursor) :: cursor

    cursor%first = 1
    cursor%last = g%ne
    cursor%current = 0
    cursor%edge = [g%edges(1)%get_entry(1), g%edges(2)%get_entry(1)]

end function coo_make_cursor



!--------------------------------------------------------------------------!
subroutine coo_get_edges(g, edges, cursor, num_edges, num_returned)        !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned

    associate(current => cursor%current)

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return; we'll either
    ! return how many edges the user asked for, or, if that amount would
    ! go beyond the final edge that the cursor is allowed to access, the
    ! all of the remaining edges
    num_returned = min(num_edges, cursor%last - cursor%current)

    ! Fill the edges array with the right slice from the graph's edges
    edges(1, 1 : num_returned) = &
                & g%edges(1)%array(current + 1 : current + num_returned)
    edges(2, 1 : num_returned) = &
                & g%edges(2)%array(current + 1 : current + num_returned)

    ! Move the cursor's current edge ahead to the last one we returned
    current = current + num_returned

    end associate

end subroutine coo_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine coo_add_edge(g, i, j)                                           !
!--------------------------------------------------------------------------!
!     NOTE: This procedure does no checking to guarantee that `i` and `j`  !
! are not already connected. This unsafe behavior is necessary in order to !
! have O(1) edge insert time for a COO graph and O(ne) storage space.      !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: i,j

    call g%edges(1)%push(i)
    call g%edges(2)%push(j)

    ! Increase the number of edges and the graph capacity
    g%ne = g%ne+1

end subroutine coo_add_edge



!--------------------------------------------------------------------------!
subroutine coo_delete_edge(g, i, j)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: i, j
    ! local variables
    integer :: k, it, jt
    type(dynamic_array) :: edges(2)

    call edges(1)%init(capacity = 4)
    call edges(2)%init(capacity = 4)

    do k = 1, g%ne
        it = g%edges(1)%get_entry(k)
        jt = g%edges(2)%get_entry(k)

        if (it /= i .or. jt /= j) then
            call edges(1)%push(it)
            call edges(2)%push(jt)
        endif
    enddo

    g%edges = edges
    g%ne = g%edges(1)%length

end subroutine coo_delete_edge



!--------------------------------------------------------------------------!
subroutine coo_graph_left_permute(g, p, edge_p)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, k

    do k = 1, g%ne
        i = g%edges(1)%get_entry(k)
        if (i /= 0) call g%edges(1)%set_entry(k, p(i))
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine coo_graph_left_permute



!--------------------------------------------------------------------------!
subroutine coo_graph_right_permute(g, p, edge_p)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, k

    do k = 1, g%ne
        i = g%edges(2)%get_entry(k)
        if (i /= 0) call g%edges(2)%set_entry(k, p(i))
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine coo_graph_right_permute



!--------------------------------------------------------------------------!
subroutine coo_destroy(g)                                                  !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(inout) :: g

    deallocate(g%edges(1)%array, g%edges(2)%array)

    g%n = 0
    g%m = 0
    g%ne = 0

end subroutine coo_destroy



!--------------------------------------------------------------------------!
subroutine coo_dump_edges(g, edges)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: k

    do k=1,g%ne
        edges(1,k) = g%edges(1)%get_entry(k)
        edges(2,k) = g%edges(2)%get_entry(k)
    enddo

end subroutine coo_dump_edges



end module coo_graphs
