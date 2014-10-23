module ll_graphs

use graph_interfaces
use types, only: dynamic_array

implicit none


!--------------------------------------------------------------------------!
type, extends(graph_interface) :: ll_graph                                 !
!--------------------------------------------------------------------------!
    ! List of `n` lists, each of which stores the neighbors of vertex `i`
    type(dynamic_array), allocatable :: lists(:)

    ! Total number of edges and maximum degree of all vertices of the graph
    integer, private :: ne, max_d
contains
    !--------------
    ! Constructors
    procedure :: init_empty => ll_graph_init
    procedure :: init_from_iterator => ll_graph_init_from_iterator

    !-----------
    ! Accessors
    procedure :: get_num_edges => ll_get_num_edges
    procedure :: get_degree => ll_get_degree
    procedure :: get_max_degree => ll_get_max_degree
    procedure :: get_neighbors => ll_get_neighbors
    procedure :: connected => ll_connected
    procedure :: find_edge => ll_find_edge
    procedure, nopass :: is_get_neighbors_fast => get_neighbors_is_fast

    !---------------
    ! Edge iterator
    procedure :: make_cursor => ll_make_cursor
    procedure :: get_edges => ll_get_edges

    !----------
    ! Mutators
    procedure :: add_edge => ll_add_edge
    procedure :: delete_edge => ll_delete_edge
    procedure :: left_permute => ll_graph_left_permute
    procedure :: right_permute => ll_graph_right_permute

    !-------------
    ! Destructors
    procedure :: destroy => ll_destroy

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => ll_dump_edges
end type ll_graph





contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_graph_init(g, n, m)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    ! local variables
    integer :: k

    call g%add_reference()

    g%n = n
    allocate(g%lists(n))

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! If we know how many neighbors each vertex has, initialize each
    ! dynamic array with the requisite amount of storage.
    do k = 1, n
        call g%lists(k)%init(capacity = 4, min_capacity = 2)
    enddo

    ! The total number of edges and max degree at initialization is zero
    g%ne = 0
    g%max_d = 0

end subroutine ll_graph_init



!--------------------------------------------------------------------------!
subroutine ll_graph_init_from_iterator(g, n, m, &                          !
                                            & iterator, get_cursor, trans) !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n, m
    procedure(edge_iterator) :: iterator
    procedure(new_cursor) :: get_cursor
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ord(2)
    integer :: i, j, k
    integer :: num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    ord = [1, 2]
    if (present(trans)) then
        if (trans) ord = [2, 1]
    endif

    call g%init(n, m)

    cursor = get_cursor()

    do while (.not. cursor%done())
        call iterator(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(ord(1), k)
            j = edges(ord(2), k)

            call g%add_edge(i, j)
        enddo
    enddo

end subroutine ll_graph_init_from_iterator




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ll_get_num_edges(g) result(k)                                     !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(in) :: g
    integer :: k

    k = g%ne

end function ll_get_num_edges



!--------------------------------------------------------------------------!
function ll_get_degree(g, i) result(d)                                     !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d

    d = g%lists(i)%length

end function ll_get_degree



!--------------------------------------------------------------------------!
function ll_get_max_degree(g) result(k)                                    !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(in) :: g
    integer :: k

    k = g%max_d

end function ll_get_max_degree



!--------------------------------------------------------------------------!
subroutine ll_get_neighbors(g, neighbors, i)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: k

    neighbors = 0
    do k = 1, g%lists(i)%length
        neighbors(k) = g%lists(i)%get_entry(k)
    enddo

end subroutine ll_get_neighbors



!--------------------------------------------------------------------------!
function ll_connected(g, i, j)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: ll_connected
    ! local variables
    integer :: k

    ll_connected = .false.
    do k = 1, g%lists(i)%length
        if (g%lists(i)%get_entry(k) == j) then
            ll_connected = .true.
            exit
        endif
    enddo

end function ll_connected



!--------------------------------------------------------------------------!
function ll_find_edge(g, i, j)                                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i, j
    integer :: ll_find_edge
    ! local variables
    integer :: k,total

    ll_find_edge = -1

    total = 0
    do k = 1, i - 1
        total = total + g%lists(k)%length
    enddo

    do k = 1, g%lists(i)%length
        if (g%lists(i)%get_entry(k) == j) then
            ll_find_edge = total + k
            exit
        endif
    enddo

end function ll_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ll_make_cursor(g) result(cursor)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: k

    cursor%first = 1
    cursor%last = g%ne
    cursor%current = 0

    k = 1
    do while (g%lists(k)%length == 0)
        k = k + 1
    enddo

    cursor%edge = [k, g%lists(k)%get_entry(1)]
    cursor%idx = 0

end function ll_make_cursor



!--------------------------------------------------------------------------!
subroutine ll_get_edges(g, edges, cursor, num_edges, num_returned)         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, k, num_added, num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return
    num_returned = min(num_edges, cursor%last - cursor%current)

    ! Find the last row we left off at
    i = cursor%edge(1)

    num_added = 0
    do while(num_added < num_returned)
        ! Either we are returning all the edges for the current vertex `i`,
        ! or the current request doesn't call for so many edges and we are
        ! returning fewer
        num_from_this_row = min(g%lists(i)%length - cursor%idx, &
                                & num_returned - num_added)

        do k = 1, num_from_this_row
            edges(1, num_added + k) = i
            edges(2, num_added + k) = g%lists(i)%get_entry(cursor%idx + k)
        enddo

        ! If we returned all nodes neighboring this vertex, increment
        ! the vertex
        !TODO replace this with bit-shifting magic
        if (num_from_this_row == g%lists(i)%length - cursor%idx) then
            i = i + 1
            cursor%idx = 0
        else
            cursor%idx = cursor%idx + num_from_this_row
        endif

        ! Increase the number of edges added
        num_added = num_added + num_from_this_row
    enddo

    cursor%current = cursor%current + num_returned
    cursor%edge(1) = i

end subroutine ll_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_add_edge(g, i, j)                                            !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i, j

    if (.not. g%connected(i, j)) then
        ! Push the new neighbor node j onto the list of i's neighbors
        call g%lists(i)%push(j)

        ! Change the graph's max degree if need be
        g%max_d = max(g%max_d, g%lists(i)%length)

        ! Increment the number of edges and graph capacity
        g%ne = g%ne + 1
    endif

end subroutine ll_add_edge



!--------------------------------------------------------------------------!
subroutine ll_delete_edge(g, i, j)                                         !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i, j
    integer :: k, jt, degree

    if (g%connected(i, j)) then
        ! Record the degree of vertex i and capacity
        degree = g%lists(i)%length

        ! Pop from the list of i's neighbors
        jt = g%lists(i)%pop()

        ! If the vertex jt popped from i's neighbors is not vertex j,
        if (jt /= j) then
            ! find where vertex j was stored and put jt there.
            do k = 1, g%lists(i)%length
                if (g%lists(i)%get_entry(k) == j) then
                    call g%lists(i)%set_entry(k, jt)
                endif
            enddo
        endif

        ! If the degree of vertex i was the max degree of the graph,
        ! check that the max degree of the graph hasn't decreased.
        !!Make this a guaranteed O(1) operation somehow
        if (degree == g%max_d) then
            g%max_d = 0
            do k=1, g%n
                g%max_d = max(g%max_d, g%lists(k)%length)
            enddo
        endif

        ! Decrement the number of edges of the graph
        g%ne = g%ne - 1
    endif

end subroutine ll_delete_edge



!--------------------------------------------------------------------------!
subroutine ll_graph_left_permute(g, p, edge_p)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, j, k
    type(dynamic_array) :: lists(g%n)

    ! If the user wants to see the edge permutation, set up the first
    ! constituent array before performing the permutation
    if (present(edge_p)) then
        allocate(edge_p(3, g%n))
        edge_p(1, 1) = 1
        do i=1, g%n-1
            edge_p(1, i + 1) = edge_p(1, i) + g%lists(i)%length
        enddo
    endif

    ! Copy all of the structure of g to a new list of lists and empty g
    do i=1, g%n
        call lists(i)%init(capacity = g%lists(i)%capacity, min_capacity = 2)
        do k = 1, g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            call lists(i)%push(j)
        enddo
        call g%lists(i)%free()
    enddo

    ! Refill g with the appropriately permuted structure
    do i = 1, g%n
        call g%lists(p(i))%init(capacity = lists(i)%capacity, min_capacity = 2)
        do k = 1, lists(i)%length
            j = lists(i)%get_entry(k)
            call g%lists(p(i))%push(j)
        enddo
        call lists(i)%free()
    enddo

    ! If the user wants to see the edge permutation, finish the remaining
    ! constituent arrays
    if (present(edge_p)) then
        edge_p(3, 1) = 1
        do i = 1, g%n-1
            edge_p(3, i + 1) = edge_p(3,i) + g%lists(i)%length
        enddo

        do i = 1, g%n
            edge_p(2, i) = edge_p(3, p(i))
        enddo

        do i = 1, g%n
            edge_p(3, i) = g%lists(p(i))%length
        enddo
    endif

end subroutine ll_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ll_graph_right_permute(g, p, edge_p)                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, j, k

    do i = 1, g%n
        do k = 1, g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            if (j /= 0) call g%lists(i)%set_entry(k, p(j))
        enddo
    enddo

    if (present(edge_p)) allocate(edge_p(0, 0))

end subroutine ll_graph_right_permute




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_destroy(g)                                                   !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer :: i

    do i = 1, g%n
        call g%lists(i)%free()
    enddo

    deallocate(g%lists)

    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_d = 0
    g%reference_count = 0

end subroutine ll_destroy




!==========================================================================!
!==== Testing, debugging & I/O                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_dump_edges(g, edges)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,k,next

    next = 0
    do i = 1, g%n
        do k = 1, g%lists(i)%length
            next = next + 1
            edges(1, next) = i
            edges(2, next) = g%lists(i)%get_entry(k)
        enddo
    enddo

end subroutine ll_dump_edges






end module ll_graphs
