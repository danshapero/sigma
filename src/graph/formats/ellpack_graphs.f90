module ellpack_graphs

use graph_interfaces
use util

implicit none


!--------------------------------------------------------------------------!
type, extends(graph_interface) :: ellpack_graph                            !
!--------------------------------------------------------------------------!
    ! Column `i` of the array `node` stores the neighbors of vertex `i`,
    ! with duplicates to fill out the remaining entries
    integer, allocatable :: node(:,:)

    ! Degrees of every vertex and max degree of the graph
    integer, allocatable :: degrees(:)
    integer :: max_d

    ! number of edges of the graph
    integer, private :: ne
contains
    !--------------
    ! Constructors
    procedure :: init_empty => ellpack_graph_init
    procedure :: init_from_iterator => ellpack_graph_init_from_iterator

    !-----------
    ! Accessors
    procedure :: get_num_edges => ellpack_get_num_edges
    procedure :: get_degree => ellpack_get_degree
    procedure :: get_max_degree => ellpack_get_max_degree
    procedure :: get_neighbors => ellpack_get_neighbors
    procedure :: connected => ellpack_connected
    procedure :: find_edge => ellpack_find_edge
    procedure, nopass :: is_get_neighbors_fast => get_neighbors_is_fast

    !---------------
    ! Edge iterator
    procedure :: make_cursor => ellpack_make_cursor
    procedure :: get_edges => ellpack_get_edges

    !----------
    ! Mutators
    procedure :: add_edge => ellpack_add_edge
    procedure :: delete_edge => ellpack_delete_edge
    procedure :: left_permute => ellpack_graph_left_permute
    procedure :: right_permute => ellpack_graph_right_permute

    !-------------
    ! Destructors
    procedure :: destroy => ellpack_destroy

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => ellpack_dump_edges

    !--------------------
    ! Auxiliary routines
    procedure, private :: add_edge_with_reallocation
end type ellpack_graph


private :: add_edge_with_reallocation



contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_graph_init(g, n, m)                                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m

    call g%add_reference()

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    allocate(g%degrees(g%n), g%node(1, g%n))
    g%degrees = 0
    g%node = 1
    g%max_d = 0
    g%ne = 0

end subroutine ellpack_graph_init



!--------------------------------------------------------------------------!
subroutine ellpack_graph_init_from_iterator(g, n, m, &                     !
                                            & iterator, get_cursor, trans) !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n, m
    procedure(edge_iterator) :: iterator
    procedure(new_cursor) :: get_cursor
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ord(2)
    integer :: i, j, k, d
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

    allocate(g%degrees(g%n))
    g%degrees = 0

    do while (.not. cursor%done())
        call iterator(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(ord(1), k)
            g%degrees(i) = g%degrees(i) + 1
        enddo
    enddo

    g%max_d = maxval(g%degrees)
    allocate(g%node(g%max_d, g%n))
    g%node = 0

    ! We're using it as a temporary for now, fear not it'll be restored
    g%degrees = 0

    cursor = get_cursor()

    do while (.not. cursor%done())
        call iterator(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(ord(1), k)
            j = edges(ord(2), k)

            if (.not. g%connected(i, j)) then
                d = g%degrees(i)
                g%node(d+1:, i) = j
                g%degrees(i) = g%degrees(i) + 1
            endif
        enddo
    enddo

end subroutine ellpack_graph_init_from_iterator




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_get_num_edges(g) result(k)                                !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer :: k

    k = g%ne

end function ellpack_get_num_edges



!--------------------------------------------------------------------------!
function ellpack_get_degree(g, i) result(d)                                !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d

    d = g%degrees(i)

end function ellpack_get_degree



!--------------------------------------------------------------------------!
function ellpack_get_max_degree(g) result(k)                               !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer :: k

    k = g%max_d

end function ellpack_get_max_degree



!--------------------------------------------------------------------------!
subroutine ellpack_get_neighbors(g, neighbors, i)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: d

    neighbors = 0
    d = g%degrees(i)

    neighbors(1:d) = g%node(1:d, i)

end subroutine ellpack_get_neighbors



!--------------------------------------------------------------------------!
function ellpack_connected(g, i, j)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i, j
    logical :: ellpack_connected
    ! local variables
    integer :: k, d

    ellpack_connected = .false.

    d = g%degrees(i)
    do k = 1, d
        if (g%node(k, i) == j) then
            ellpack_connected = .true.
            exit
        endif
    enddo

end function ellpack_connected



!--------------------------------------------------------------------------!
function ellpack_find_edge(g, i, j)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i, j
    integer :: ellpack_find_edge
    ! local variables
    integer :: k, d

    ellpack_find_edge = -1

    d = g%degrees(i)
    do k = 1, d
        if (g%node(k, i) == j) ellpack_find_edge = g%max_d * (i - 1) + k
    enddo

end function ellpack_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_make_cursor(g) result(cursor)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: i

    cursor%first = 1
    cursor%last = g%ne
    cursor%current = 0

    i = 1
    do while(g%degrees(i) == 0)
        i = i + 1
    enddo
    cursor%edge = [i, g%node(1, i)]
    cursor%idx = 0

end function ellpack_make_cursor



!--------------------------------------------------------------------------!
subroutine ellpack_get_edges(g, edges, cursor, num_edges, num_returned)    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: num_edges
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, k, num, bt

    associate(idx => cursor%idx)

    ! Set up the returned edges to be 0
    edges = 0

    ! Find out how many nodes' edges the current request encompasses
    num_returned = min(num_edges, cursor%last - cursor%current)

    ! Find the last vertex that the edge iterator left off at
    i = cursor%edge(1)

    ! Total number of edges from the current batch that we've fetched so far
    k = 0

    ! Keep fetching more edges as long as we haven't got the whole batch
    do while(k < num_returned)
        ! Number of edges we'll return that neighbor vertex `i`.
        ! Either we're returning the rest of the row, i.e.
        !     g%degrees(i) - idx
        ! edges, or the user hasn't asked for that many edges, so they get
        !     num_returned - k.
        num = min(g%degrees(i) - idx, num_returned - k)

        edges(1, k+1 : k+num) = i
        edges(2, k+1 : k+num) = g%node(idx+1 : idx+num, i)

        ! The following statements are some bit-shifting magic in order to
        ! avoid the conditional
        !     if (num == g%degrees(i) - idx) then
        !         i = i + 1
        !         idx = 0
        !     else
        !         idx = idx + num
        !     endif
        ! The two are completely equivalent, but we'd like to avoid
        ! branching where possible.
        bt = (sign(1, num - (g%degrees(i) - idx)) + 1) / 2
        i = i + bt
        idx = (1 - bt) * (idx + num)

        k = k + num
    enddo

    cursor%current = cursor%current + num_returned
    cursor%edge(1) = i

    end associate

end subroutine ellpack_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_add_edge(g, i, j)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i, j
    ! local variables
    integer :: k

    ! If vertices i,j are already connected, we needn't add the edge
    if (.not.g%connected(i, j)) then
        ! Find the index in g%node(:,i) where we can add in node j
        k = g%degrees(i)

        ! If there is room to add j, then do so
        if (k < g%max_d) then
            ! Set the entire rest of the row in the array `node` to be
            ! equal to `j`. That way, the edge iterator never returns `0`,
            ! only duplicates of edges that already exist.
            g%node(k + 1 :, i) = j

            ! Increment the degree count for vertex `i`
            g%degrees(i) = g%degrees(i) + 1

            ! Increment the total number of edges of the graph
            g%ne = g%ne + 1

        ! If there is no room, then degree(i) = max degree of g and we need
        ! to reallocate space. This is deferred to a helper subroutine at
        ! the end of this module.
        else
            call g%add_edge_with_reallocation(i, j)
        endif
    endif

end subroutine ellpack_add_edge



!--------------------------------------------------------------------------!
subroutine ellpack_delete_edge(g, i, j)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: indx, d
    integer, allocatable :: node(:, :)
    logical :: max_degree_decrease

    ! If nodes i,j are not connected to begin with, there is no edge to
    ! delete and thus nothing to do
    if (g%connected(i,j)) then
        ! Record what the max degree of `g` was before removing the edge
        d = g%max_d

        ! Set a boolean to be true if we're removing an edge from a node
        ! of maximum degree
        max_degree_decrease = (g%degrees(i) == d)

        ! Find the location indx in memory where edge (i,j) is stored
        do indx = 1, g%max_d
            if (g%node(indx, i) == j) exit
        enddo

        ! Overwrite indx with the other nodes connected to i
        g%node(indx : g%max_d - 1, i) = g%node(indx + 1 : g%max_d, i)

        ! Replace the last vertex connected to `i` with a copy of the
        ! previous edge. This is so there are no `0` entries in the array
        ! `node`.
        if (g%max_d > 0) then
            g%node(g%max_d, i) = g%node(g%max_d - 1, i)
        endif

        ! Decrement the number of edges
        g%ne = g%ne - 1

        ! Decrement the degree of node `i`
        g%degrees(i) = g%degrees(i) - 1

        ! If node i had max degree, check all the other nodes to see if
        ! the max degree has decreased
        if (max_degree_decrease) then
            g%max_d = maxval(g%degrees)

            ! If the max degree really has decreased, we need to reduce
            ! the storage space for `g`
            if (g%max_d < d) then
                ! Make a temporary array of size sufficient for the reduced
                ! number of edges of `g`
                allocate(node(g%max_d, g%n))

                ! Copy the parts of `g` that still remain to the temporary
                node(1 : g%max_d, :) = g%node(1 : g%max_d, :)

                ! Move the allocation status from the temporary to `g`
                call move_alloc(from = node, to = g%node)
            endif
        endif

    endif

end subroutine ellpack_delete_edge



!--------------------------------------------------------------------------!
subroutine ellpack_graph_left_permute(g, p, edge_p)                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, node(g%max_d, g%n)

    ! Copy the list of neighbors of vertex `i` from `g` to the temporary
    ! array `node` in column `p(i)`
    do i = 1, g%n
        node(:, p(i)) = g%node(:, i)
    enddo

    ! Copy the temporary array into `g`
    g%node = node

    ! Update the degrees of all the vertices
    g%degrees(p) = g%degrees

    if (present(edge_p)) then
        ! Report the resulting edge permutation
        allocate(edge_p(3, g%n))

        do i=1,g%n
            edge_p(1, i) = g%max_d * (i - 1) + 1
            edge_p(2, i) = g%max_d * (p(i) - 1) + 1
            edge_p(3, i) = g%max_d
        enddo
    endif

end subroutine ellpack_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_graph_right_permute(g, p, edge_p)                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i, j, k

    do i = 1, g%n
        do k = 1, g%max_d
            j = g%node(k, i)
            if (j /= 0) g%node(k, i) = p(j)
        enddo
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine ellpack_graph_right_permute




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_destroy(g)                                              !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g

    deallocate(g%node, g%degrees)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_d = 0
    g%reference_count = 0

end subroutine ellpack_destroy




!==========================================================================!
!==== Testing, debugging & I/O                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_dump_edges(g, edges)                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i, j, k, d, next

    next = 0
    do i = 1, g%n
        d = g%degrees(i)
        do k = 1, d
            j = g%node(k, i)
            if (j /= 0) then
                next = next + 1
                edges(1, next) = i
                edges(2, next) = j
            else
                exit
            endif
        enddo
    enddo

end subroutine ellpack_dump_edges



!--------------------------------------------------------------------------!
subroutine add_edge_with_reallocation(g, i, j)                             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i, j
    ! local variables
    integer :: k
    integer, allocatable :: node(:, :)

    ! Make a temporary array sufficient for the new size of `g`
    allocate(node(g%max_d + 1, g%n))
    node = 0

    ! Copy all the old entries from `g` into the temporary array
    node(1 : g%max_d, :) = g%node(1 : g%max_d, :)

    ! Get rid of all `0` entries in the array by putting in copies of 
    ! already-existing edges
    if (g%max_d > 0) then
        do k = 1, g%n
            node(g%max_d + 1, k) = node(g%max_d, k)
        enddo
    endif

    ! Put vertex `j` in the list of neighbors of vertex `i`
    node(g%max_d + 1, i) = j

    ! Move the allocation status from the temporary to `g`
    call move_alloc(from = node, to = g%node)

    ! Increment the number of edges of `g`
    g%ne = g%ne + 1

    ! Increment the degree of vertex `i`
    g%degrees(i) = g%degrees(i) + 1

    ! Increment the max degree of `g`
    g%max_d = g%max_d + 1

end subroutine add_edge_with_reallocation





end module ellpack_graphs

