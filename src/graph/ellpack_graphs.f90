module ellpack_graphs

use graph_interface
use util

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ellpack_graph                                      !
!--------------------------------------------------------------------------!
    integer, allocatable :: node(:,:), degrees(:)
contains
    !--------------
    ! Constructors
    procedure :: init => ellpack_graph_init
    procedure :: copy => ellpack_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => ellpack_degree
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

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    allocate(g%degrees(g%n), g%node(1, g%n))
    g%degrees = 0
    g%node = 1
    g%max_degree = 0
    g%ne = 0

end subroutine ellpack_graph_init



!--------------------------------------------------------------------------!
subroutine ellpack_graph_copy(g, h, trans)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    class(graph), intent(in)            :: h
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ind(2), order(2), nv(2), k
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    nv = [h%n, h%m]
    order = [1, 2]

    ! Check if we're copying h or h with all directed edges reversed
    if (present(trans)) then
        if (trans) then
            nv = [h%m, h%n]
            order = [2, 1]
        endif
    endif

    ! Copy all the attributes of g from those of h
    g%n = nv(1)
    g%m = nv(2)
    g%ne = 0
    g%max_degree = h%max_degree

    ! Allocate space for the main node array of g
    allocate(g%node(g%max_degree, g%n), g%degrees(g%n))
    g%node = 0

    ! Make an edge iterator for the copied graph h
    cursor = h%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    ! Iterate through all the edges of h
    do n = 1, num_batches
        ! Get a chunk of edges from h
        call h%get_edges(edges, cursor, batch_size, num_returned)

        ! Add each edge from the chunk into g
        do k = 1, num_returned
            ind = edges(order, k)
            call g%add_edge(ind(1),ind(2))
            g%degrees(ind(1)) = g%degrees(ind(1)) + 1
        enddo
    enddo

end subroutine ellpack_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_degree(g, i) result(d)                                    !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d

    d = g%degrees(i)

end function ellpack_degree



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
    integer :: k

    ellpack_connected = .false.

    do k = 1, g%max_degree
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
    integer :: k

    ellpack_find_edge = -1

    do k = g%max_degree, 1, -1
        if (g%node(k, i) == j) then
            ellpack_find_edge = (i - 1) * g%max_degree + k
            exit
        endif
    enddo

end function ellpack_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_make_cursor(g) result(cursor)                             !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0
    cursor%edge = [1, g%node(1, 1)]
    cursor%indx = 0

end function ellpack_make_cursor



!--------------------------------------------------------------------------!
subroutine ellpack_get_edges(g, edges, cursor, num_edges, num_returned)    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, d, i1, i2, num_added, num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Find out how many nodes' edges the current request encompasses
    num_returned = min(num_edges, cursor%final - cursor%current)

    ! Find the starting and ending nodes for this edge retrieval
    i1 = cursor%current / g%max_degree + 1
    i2 = (cursor%current + num_returned - 1) / g%max_degree + 1

    ! Set the number of edges added to 0
    num_added = 0

    ! Loop from the starting node to the ending node
    do i = i1, i2
        ! Find how many edges we're retrieving from this row
        num_from_this_row = min(g%max_degree - cursor%indx, &
                                & num_returned - num_added)

        ! Fill in the return array
        edges(1, num_added + 1 : num_added + num_from_this_row) = i
        edges(2, num_added + 1 : num_added + num_from_this_row) = &
            & g%node(cursor%indx + 1 : cursor%indx + num_from_this_row, i)

        ! Increment the number of edges added
        num_added = num_added + num_from_this_row

        ! Modify the index storing the place within the row that we left off
        cursor%indx = mod(cursor%indx + num_from_this_row, g%max_degree)
    enddo

    cursor%current = cursor%current + num_returned

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
        if (k < g%max_degree) then
            ! Set the entire rest of the row in the array `node` to be
            ! equal to `j`. That way, the edge iterator never returns `0`,
            ! only duplicates of edges that already exist.
            g%node(k + 1 : g%max_degree, i) = j

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
        d = g%max_degree

        ! Set a boolean to be true if we're removing an edge from a node
        ! of maximum degree
        max_degree_decrease = (g%degrees(i) == d)

        ! Find the location indx in memory where edge (i,j) is stored
        do indx = 1, g%max_degree
            if (g%node(indx, i) == j) exit
        enddo

        ! Overwrite indx with the other nodes connected to i
        g%node(indx : g%max_degree - 1, i) &
                                    & = g%node(indx + 1 : g%max_degree, i)

        ! Replace the last vertex connected to `i` with a copy of the
        ! previous edge. This is so there are no `0` entries in the array
        ! `node`.
        if (g%max_degree > 0) then
            g%node(g%max_degree, i) = g%node(g%max_degree - 1, i)
        endif

        ! Decrement the number of edges
        g%ne = g%ne - 1

        ! Decrement the degree of node `i`
        g%degrees(i) = g%degrees(i) - 1

        ! If node i had max degree, check all the other nodes to see if
        ! the max degree has decreased
        if (max_degree_decrease) then
            g%max_degree = maxval(g%degrees)

            ! If the max degree really has decreased, we need to reduce
            ! the storage space for `g`
            if (g%max_degree < d) then
                ! Make a temporary array of size sufficient for the reduced
                ! number of edges of `g`
                allocate(node(g%max_degree, g%n))

                ! Copy the parts of `g` that still remain to the temporary
                node(1 : g%max_degree, :) = g%node(1 : g%max_degree, :)

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
    integer :: i, node(g%max_degree, g%n)

    do i=1,g%n
        node(:, p(i)) = g%node(:, i)
    enddo

    g%node = node

    if (present(edge_p)) then
        ! Report the resulting edge permutation
        allocate(edge_p(3, g%n))

        do i=1,g%n
            edge_p(1, i) = g%max_degree * (i - 1) + 1
            edge_p(2, i) = g%max_degree * (p(i) - 1) + 1
            edge_p(3, i) = g%max_degree
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
        do k = 1, g%max_degree
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

    deallocate(g%node)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0

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
    allocate(node(g%max_degree + 1, g%n))
    node = 0

    ! Copy all the old entries from `g` into the temporary array
    node(1 : g%max_degree, :) = g%node(1 : g%max_degree, :)

    ! Get rid of all `0` entries in the array by putting in copies of 
    ! already-existing edges
    do k = 1, g%n
        g%node(g%max_degree + 1, k) = g%node(g%max_degree, k)
    enddo

    ! Put vertex `j` in the list of neighbors of vertex `i`
    node(g%max_degree + 1, i) = j

    ! Move the allocation status from the temporary to `g`
    call move_alloc(from = node, to = g%node)

    ! Increment the number of edges of `g`
    g%ne = g%ne + 1

    ! Increment the degree of vertex `i`
    g%degrees(i) = g%degrees(i) + 1

    ! Increment the max degree of `g`
    g%max_degree = g%max_degree + 1

end subroutine add_edge_with_reallocation





end module ellpack_graphs

