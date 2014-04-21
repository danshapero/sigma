module ellpack_graphs

use graphs
use util

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ellpack_graph                                      !
!--------------------------------------------------------------------------!
    integer, allocatable :: node(:,:)
    integer :: max_neighbors
contains
    !--------------
    ! Constructors
    procedure :: init_const_degree => ellpack_init_const_degree
    procedure :: init_variable_degree => ellpack_init_variable_degree
    procedure :: copy => ellpack_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => ellpack_degree
    procedure :: neighbors => ellpack_neighbors
    procedure :: connected => ellpack_connected
    procedure :: find_edge => ellpack_find_edge

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
    procedure :: compress => ellpack_graph_compress
    procedure :: decompress => ellpack_graph_decompress

    !-------------
    ! Destructors
    procedure :: free => ellpack_free

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => ellpack_dump_edges

    !--------------------
    ! Auxiliary routines
    procedure, private :: add_edge_with_reallocation
    procedure :: ellpack_max_degree_decrease
end type ellpack_graph





contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_init_const_degree(g,n,m,degree)                         !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, degree

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    if (present(degree)) then
        g%max_neighbors = degree
    else
        g%max_neighbors = 1
    endif

    allocate(g%node(g%max_neighbors,g%n))
    g%node = 0
    g%max_degree = 0
    g%ne = 0
    g%capacity = g%max_neighbors*g%n

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine ellpack_init_const_degree



!--------------------------------------------------------------------------!
subroutine ellpack_init_variable_degree(g,n,m,degrees)                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: n, degrees(:)
    integer, intent(in), optional :: m

    g%n = n

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    g%max_neighbors = maxval(degrees)

    allocate(g%node(g%max_neighbors,g%n))
    g%node = 0
    g%max_degree = 0
    g%ne = 0
    g%capacity = g%max_neighbors*g%n

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine ellpack_init_variable_degree



!--------------------------------------------------------------------------!
subroutine ellpack_graph_copy(g,h,trans)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    class(graph), intent(in)            :: h
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ind(2),order(2),nv(2),k,n,num_blocks,num_returned,edges(2,64)
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

    ! Mark the graph as mutable
    g%mutable = .true.

    ! Copy all the attributes of g from those of h
    g%n = nv(1)
    g%m = nv(2)
    g%ne = 0
    g%max_degree = h%max_degree
    g%max_neighbors = g%max_degree
    g%capacity = g%max_degree * g%n

    ! Allocate space for the main node array of g
    allocate(g%node(g%max_degree,g%n))
    g%node = 0

    ! Make an edge iterator for the copied graph h
    cursor = h%make_cursor(0)
    num_blocks = (cursor%final-cursor%current)/64+1

    ! Iterate through all the edges of h
    do n=1,num_blocks
        ! Get a chunk of edges from h
        edges = h%get_edges(cursor,64,num_returned)

        ! Add each edge from the chunk into g
        do k=1,num_returned
            ind = edges(order,k)

            if (ind(1)/=0 .and. ind(2)/=0) then
                call g%add_edge(ind(1),ind(2))
            endif
        enddo
    enddo

end subroutine ellpack_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_degree(g,i) result(d)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d
    ! local variables
    integer :: j, k

    d = 0
    do k=1,g%max_degree
        j = g%node(k,i)
        if (j/=0) d = k
    enddo

end function ellpack_degree



!--------------------------------------------------------------------------!
subroutine ellpack_neighbors(g,i,nbrs)                                     !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i
    integer, intent(out) :: nbrs(:)

    nbrs = 0
    nbrs(1:g%max_degree) = g%node(:,i)

end subroutine ellpack_neighbors



!--------------------------------------------------------------------------!
function ellpack_connected(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: ellpack_connected
    ! local variables
    integer :: k

    ellpack_connected = .false.

    do k=1,g%max_degree
        if (g%node(k,i)==j) ellpack_connected = .true.
    enddo

end function ellpack_connected



!--------------------------------------------------------------------------!
function ellpack_find_edge(g,i,j)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: ellpack_find_edge
    ! local variables
    integer :: k

    ellpack_find_edge = -1

    do k=g%max_degree,1,-1
        if (g%node(k,i)==j) ellpack_find_edge = (i-1)*g%max_degree+k
    enddo

end function ellpack_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ellpack_make_cursor(g,thread) result(cursor)                      !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor

    cursor%start = 1
    cursor%final = g%capacity
    cursor%current = 0
    cursor%edge = [1, g%node(1,1)]
    cursor%indx = 0

end function ellpack_make_cursor



!--------------------------------------------------------------------------!
function ellpack_get_edges(g,cursor,num_edges,num_returned) result(edges)  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    integer :: edges(2,num_edges)
    ! local variables
    integer :: i,i1,i2,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Find out how many nodes' edges the current request encompasses
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! Find the starting and ending nodes for this edge retrieval
    i1 = cursor%current/g%max_neighbors+1
    i2 = (cursor%current+num_returned-1)/g%max_neighbors+1

    ! Set the number of edges added to 0
    num_added = 0

    ! Loop from the starting node to the ending node
    do i=i1,i2
        ! Find how many edges we're retrieving from this row
        num_from_this_row = min(g%max_degree-cursor%indx, &
                                & num_returned-num_added)

        ! Fill in the return array
        edges(1,num_added+1:num_added+num_from_this_row) = i
        edges(2,num_added+1:num_added+num_from_this_row) = &
            & g%node(cursor%indx+1:cursor%indx+num_from_this_row,i)

        ! Increment the number of edges added
        num_added = num_added+num_from_this_row

        ! Modify the index storing the place within the row that we left off
        cursor%indx = mod(cursor%indx+num_from_this_row,g%max_degree)
    enddo

    cursor%current = cursor%current+num_returned

end function ellpack_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_add_edge(g,i,j)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k,indx

    if (.not.g%mutable) then
        print *, 'Attempting to add edge to an immutable ellpack graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! If vertices i,j are already connected, we needn't add the edge
    if (.not.g%connected(i,j)) then
        ! Find the index in g%node(:,i) where we can add in node j
        indx = -1

        do k=g%max_neighbors,1,-1
            if (g%node(k,i)==0) then
                indx = k
            endif
        enddo

        ! If there is room to add j, then do so
        if (indx/=-1) then
            g%node(indx,i) = j
            g%ne = g%ne+1
            g%max_degree = max(g%max_degree,indx)
        ! If there is no room, then degree(i) = max degree of g.
        else
            call g%add_edge_with_reallocation(i,j)
        endif
    endif

end subroutine ellpack_add_edge



!--------------------------------------------------------------------------!
subroutine ellpack_delete_edge(g,i,j)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: indx
    logical :: max_degree_decrease

    if (.not.g%mutable) then
        print *, 'Attempted to delete edge from immutable ellpack graph'
        print *, 'Terminating.'
        call exit(1)
    else
        ! If nodes i,j are not connected to begin with, there is no edge to
        ! delete and thus nothing to do
        if (g%connected(i,j)) then
            ! Set a boolean to be true if we're removing an edge from a node
            ! of maximum degree
            max_degree_decrease = (g%node(g%max_degree,i)/=0)

            ! Find the location indx in memory where edge (i,j) is stored
            do indx=1,g%max_degree
                if (g%node(indx,i)==j) exit
            enddo

            ! Overwrite indx with the other nodes connected to i
            g%node(indx:g%max_degree-1,i) = g%node(indx+1:g%max_degree,i)

            ! Zero out the last node connected to i
            g%node(g%max_degree,i) = 0

            ! Decrement the number of edges
            g%ne = g%ne-1

            ! If node i had max degree, check all the other nodes to see if
            ! the max degree has decreased
            if (max_degree_decrease) then
                max_degree_decrease = (maxval(g%node(g%max_degree,:))==0)
            endif

            ! If so, decrement the max_degree member of g
            if (max_degree_decrease) g%max_degree = g%max_degree-1
        endif
    endif

end subroutine ellpack_delete_edge



!--------------------------------------------------------------------------!
subroutine ellpack_graph_left_permute(g,p,edge_p)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,node(g%max_neighbors,g%n)

    do i=1,g%n
        node(:,p(i)) = g%node(:,i)
    enddo

    g%node = node

    if (present(edge_p)) then
        ! Report the resulting edge permutation
        allocate(edge_p(3,g%n))

        do i=1,g%n
            edge_p(1,i) = g%max_neighbors*(i-1)+1
            edge_p(2,i) = g%max_neighbors*(p(i)-1)+1
            edge_p(3,i) = g%max_neighbors
        enddo
    endif

end subroutine ellpack_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ellpack_graph_right_permute(g,p,edge_p)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,j,k

    do i=1,g%n
        do k=1,g%max_degree
            j = g%node(k,i)
            if (j/=0) g%node(k,i) = p(j)
        enddo
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine ellpack_graph_right_permute



!--------------------------------------------------------------------------!
subroutine ellpack_graph_compress(g,edge_p)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, allocatable, intent(inout), optional :: edge_p(:,:)
    ! local variables
    integer :: i,j,k,jt
    integer, allocatable :: node(:,:)

    ! If the maximum number of neighbors possible for each vertex of the
    ! graph is greater than the maximum degree of the graph, we can reduce
    ! the storage needed
    if (g%max_neighbors/=g%max_degree) then
        allocate(node(g%max_degree,g%n))
        g%max_neighbors = g%max_degree
        g%capacity = g%max_degree * g%n

        ! Copy the array g%node into a smaller temporary array node
        do i=1,g%n
            node(1:g%max_degree,i) = g%node(1:g%max_degree,i)
        enddo

        ! Transfer the allocation status of node to g%node
        call move_alloc(from=node, to=g%node)

        ! Compute the edge permutation if the caller wants it
        if (present(edge_p)) then
            allocate(edge_p(3,g%n))

            do i=1,g%n
                edge_p(1,i) = g%max_neighbors*(i-1)+1
                edge_p(2,i) = g%max_degree*(i-1)+1
                edge_p(3,i) = g%max_degree
            enddo
        endif
    else
        if (present(edge_p)) allocate(edge_p(0,0))
    endif

    ! Fill out any remaining null edges with copies of existing edges so
    ! that the edge iterator never returns a null edge
    do i=1,g%n
        jt = g%node(1,i)
        do k=2,g%max_degree
            j = g%node(k,i)
            if (j/=0) then
                jt = j
            else
                g%node(k,i) = jt
            endif
        enddo
    enddo

    ! Mark the graph as immutable
    g%mutable = .false.

end subroutine ellpack_graph_compress



!--------------------------------------------------------------------------!
subroutine ellpack_graph_decompress(g)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    ! local variables
    integer :: i,j1,j2,k

    g%mutable = .true.

    ! If there were any duplicate edges created as a result of compression,
    ! replace them with null edges.
    do i=1,g%n
        do k=g%max_neighbors,2,-1
            j1 = g%node(k,i)
            j2 = g%node(k-1,i)
            if (j1==j2) g%node(k,i) = 0
        enddo
    enddo
    ! Note that if you had an isolated vertex with no edges at all, this
    ! would spuriously leave it with one edge left. But that's a pretty
    ! weird edge case.

end subroutine ellpack_graph_decompress




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_free(g)                                                 !
!--------------------------------------------------------------------------!
    class(ellpack_graph), intent(inout) :: g

    deallocate(g%node)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0
    g%capacity = 0

    g%mutable = .true.

end subroutine ellpack_free




!==========================================================================!
!==== Testing, debugging & I/O                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ellpack_dump_edges(g,edges)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,j,k,next

    next = 0
    do i=1,g%n
        do k=1,g%max_degree
            j = g%node(k,i)
            if (j/=0) then
                next = next+1
                edges(1,next) = i
                edges(2,next) = j
            else
                exit
            endif
        enddo
    enddo

end subroutine ellpack_dump_edges



!--------------------------------------------------------------------------!
subroutine add_edge_with_reallocation(g,i,j)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer, allocatable :: node(:,:)

    allocate(node(g%max_neighbors+1,g%n))
    node = 0
    node(1:g%max_neighbors,:) = g%node(1:g%max_neighbors,:)
    node(g%max_neighbors+1,i) = j

    call move_alloc(from=node, to=g%node)

    g%ne = g%ne+1
    g%max_degree = g%max_degree+1
    g%max_neighbors = g%max_neighbors+1
    g%capacity = g%n*g%max_neighbors

end subroutine add_edge_with_reallocation



!--------------------------------------------------------------------------!
subroutine ellpack_max_degree_decrease(g)                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ellpack_graph), intent(inout) :: g
    ! local variables
    integer :: k, node(g%max_degree,g%n)

    node = g%node
    deallocate(g%node)
    allocate(g%node(g%max_degree-1,g%n))

    g%node = 0
    do k=1,g%n
        g%node(1:g%max_degree-1,k) = node(1:g%max_degree-1,k)
    enddo

    g%max_degree = g%max_degree-1

end subroutine ellpack_max_degree_decrease





end module ellpack_graphs
