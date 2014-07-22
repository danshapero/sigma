module coo_graphs

use graph_interface
use types, only: dynamic_array

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: coo_graph                                          !
!--------------------------------------------------------------------------!
    type(dynamic_array) :: edges(2)
contains
    !--------------
    ! Constructors
    procedure :: init => coo_graph_init
    procedure :: copy => coo_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => coo_degree
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
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m
    ! local variables
    integer :: ne

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
    g%max_degree = 0

end subroutine coo_graph_init



!--------------------------------------------------------------------------!
subroutine coo_graph_copy(g, h, trans)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(coo_graph), intent(inout) :: g
    class(graph), intent(in)        :: h
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ind(2), ord(2), nv(2), k
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    nv = [h%n, h%m]
    ord = [1, 2]

    ! Check if we're copying h, or h with all edges reversed
    if (present(trans)) then
        if (trans) then
            nv = [h%m, h%n]
            ord = [2, 1]
        endif
    endif

    ! Copy all the attributes of h to g
    g%n = nv(1)
    g%m = nv(2)
    g%ne = h%ne
    g%max_degree = h%max_degree

    ! Allocate space in the two dynamic arrays for the edges of g
    call g%edges(1)%init(capacity = g%ne + 16)
    call g%edges(2)%init(capacity = g%ne + 16)

    ! Make an edge iterator for the copied graph h
    cursor = h%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    ! Iterate through all the edges of h
    do n = 1, num_batches
        ! Get a chunk of edges of h
        call h%get_edges(edges, cursor, batch_size, num_returned)

        ! Add each edge from the chunk into g
        do k = 1, num_returned
            ind = edges(ord, k)

            call g%edges(1)%push(ind(1))
            call g%edges(2)%push(ind(2))
        enddo
    enddo

end subroutine coo_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function coo_degree(g, i) result(d)                                        !
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

end function coo_degree



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
    integer :: j, k, l, next
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

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0
    cursor%edge = [g%edges(1)%get_entry(1), g%edges(2)%get_entry(1)]

end function coo_make_cursor



!--------------------------------------------------------------------------!
subroutine coo_get_edges(g, edges, cursor, num_edges, num_returned)        !
!--------------------------------------------------------------------------!
    class(coo_graph), intent(in) :: g
    integer, intent(out) :: edges(2, num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned

    associate(current => cursor%current)

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return; we'll either
    ! return how many edges the user asked for, or, if that amount would
    ! go beyond the final edge that the cursor is allowed to access, the
    ! all of the remaining edges
    num_returned = min(num_edges, cursor%final - cursor%current)

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
    integer :: k, it, jt, ir, jr

    ! Find an edge `(ir, jr)` distinct from `(i, j)`
    do k = 1, g%ne
        it = g%edges(1)%get_entry(k)
        jt = g%edges(2)%get_entry(k)

        if (it /= i .or. jt /= j) then
            ir = it
            jr = jt
        endif
    enddo

    ! Find any instances of the edge `(i, j)` and replace them with a
    ! duplicate of the edge `(ir, jr)`. Note that COO graphs are allowed to
    ! store duplicate edges.
    do k = 1, g%ne
        it = g%edges(1)%get_entry(k)
        jt = g%edges(2)%get_entry(k)

        if (it == i .and. jt == j) then
            call g%edges(1)%set_entry(k, ir)
            call g%edges(2)%set_entry(k, ir)
        endif
    enddo

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
    g%max_degree = 0

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
