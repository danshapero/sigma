module cs_graphs

use graph_interface
use util

implicit none



!--------------------------------------------------------------------------!
type, extends(graph) :: cs_graph                                           !
!--------------------------------------------------------------------------!
    integer, allocatable :: ptr(:), node(:)
contains
    !--------------
    ! Constructors
    procedure :: init_const_degree => cs_init_const_degree
    procedure :: init_variable_degree => cs_init_variable_degree
    procedure :: copy => cs_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => cs_degree
    procedure :: get_neighbors => cs_get_neighbors
    procedure :: connected => cs_connected
    procedure :: find_edge => cs_find_edge

    !---------------
    ! Edge iterator
    procedure :: make_cursor => cs_make_cursor
    procedure :: get_edges => cs_get_edges

    !----------
    ! Mutators
    procedure :: add_edge  => cs_add_edge
    procedure :: delete_edge => cs_delete_edge
    procedure :: left_permute => cs_graph_left_permute
    procedure :: right_permute => cs_graph_right_permute
    procedure :: compress => cs_graph_compress
    procedure :: decompress => cs_graph_decompress

    !-------------
    ! Destructors
    procedure :: free => cs_free

    !--------------------------
    ! Testing, debugging & I/O
    procedure :: dump_edges => cs_dump_edges

    !--------------------
    ! Auxiliary routines
    procedure, private :: add_edge_with_reallocation
    procedure :: sort_node
    procedure, private :: max_degree_update
end type cs_graph




contains




!==========================================================================!
!==== Constructors                                                     ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_init_const_degree(g,n,m,degree)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, degree
    ! local variables
    integer :: k, ne

    ! Set the number of (left-)nodes in the graph and allocate the node
    ! pointer array ptr
    g%n = n
    allocate( g%ptr(n+1) )

    ! If the graph is not square, set the number of right-nodes
    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! If the degree argument is present, we have an upper bound on the
    ! number of edges to allocate in the graph
    if (present(degree)) then
        ne = degree*g%n
    else
        ne = g%n
    endif

    ! Allocate the node array to be the number of edges, whether we actually
    ! have an upper bound on what that is or we've set it to some default
    allocate( g%node(ne) )
    g%node = 0
    g%capacity = ne

    ! If we know how many neighbors each vertex has, fill in the ptr
    ! array accordingly
    if (present(degree)) then
        do k=1,g%n+1
            g%ptr(k) = (k-1)*degree+1
        enddo
    else
        g%ptr = 1
    endif

    ! The total number of edges and max degree at initialization is zero
    g%ne = 0
    g%max_degree = 0

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine cs_init_const_degree



!--------------------------------------------------------------------------!
subroutine cs_init_variable_degree(g,n,m,degrees)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: n, degrees(:)
    integer, intent(in), optional :: m
    ! local variables
    integer :: k,ne

    ! Set the number of (left-)nodes in the graph and allocate the node
    ! pointer array ptr
    g%n = n
    allocate( g%ptr(n+1) )

    ! If the graph is not square, set the number of right-nodes
    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! Set the number of edges to allocate in the node array
    ne = sum(degrees)

    ! Allocate the node array to be the number of edges
    allocate( g%node(ne) )
    g%node = 0
    g%capacity = ne

    ! If we know how many neighbors each vertex has, fill in the ptr
    ! array accordingly
    g%ptr(1) = 1
    do k=1,g%n
        g%ptr(k+1) = g%ptr(k)+degrees(k)
    enddo

    ! The total number of edges and max degree at initialization is zero
    g%ne = 0
    g%max_degree = 0

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine cs_init_variable_degree



!--------------------------------------------------------------------------!
subroutine cs_graph_copy(g,h,trans)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    class(graph), intent(in)       :: h
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ind(2),order(2),nv(2),k,n,num_blocks,num_returned,edges(2,64)
    type(graph_edge_cursor) :: cursor

    nv = [h%n, h%m]
    order = [1, 2]

    ! Check whether we're copying h or h with all directed edes reversed
    if (present(trans)) then
        if (trans) then
            nv = [h%m, h%n]
            order = [2, 1]
        endif
    endif

    ! Mark the graph as mutable
    g%mutable = .true.

    ! Copy all of h's attributes to g
    g%n = nv(1)
    g%m = nv(2)
    g%ne = 0
    g%capacity = h%ne
    g%max_degree = h%max_degree

    ! Allocate g's ptr and node arrays
    allocate(g%ptr(g%n+1),g%node(g%capacity))

    ! Get a cursor from h with which to iterate through its edges
    cursor = h%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    ! Fill out the ptr array
    g%ptr = 0

    ! Iterate through the edges of h first to fill out the ptr array of g
    do n=1,num_blocks
        ! Get a chunk of edges from h
        call h%get_edges(edges,cursor,64,num_returned)

        ! For each edge,
        do k=1,num_returned
            ind = edges(order,k)

            ! If that edge isn't null
            if (ind(1)/=0 .and. ind(2)/=0) then
                ! Increment a counter
                g%ptr(ind(1)+1) = g%ptr(ind(1)+1)+1
            endif
        enddo
    enddo

    g%ptr(1) = 1
    do k=1,g%n
        g%ptr(k+1) = g%ptr(k)+g%ptr(k+1)
    enddo

    ! Iterate through the edges of h again to fill the node array of g
    cursor = h%make_cursor(0)

    g%node = 0

    do n=1,num_blocks
        call h%get_edges(edges,cursor,64,num_returned)

        do k=1,num_returned
            ind = edges(order,k)

            if (ind(1)/=0 .and. ind(2)/=0) then
                call g%add_edge(ind(1),ind(2))
            endif
        enddo
    enddo

end subroutine cs_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_degree(g,i) result(d)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d
    ! local variables
    integer :: j,k

    d = 0
    do k=g%ptr(i),g%ptr(i+1)-1
        j = g%node(k)
        if (j/=0) d = k-g%ptr(i)+1
    enddo

end function cs_degree



!--------------------------------------------------------------------------!
subroutine cs_get_neighbors(g,neighbors,i)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: start,finish,degree

    degree = g%ptr(i+1)-g%ptr(i)
    start = g%ptr(i)
    finish = g%ptr(i+1)-1

    neighbors = 0
    neighbors(1:degree) = g%node(start:finish)

end subroutine cs_get_neighbors



!--------------------------------------------------------------------------!
function cs_connected(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: cs_connected
    ! local variables
    integer :: k

    cs_connected = .false.

    do k=g%ptr(i),g%ptr(i+1)-1
        if (g%node(k)==j) cs_connected = .true.
    enddo

end function cs_connected



!--------------------------------------------------------------------------!
function cs_find_edge(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: cs_find_edge
    ! local variables
    integer :: k

    cs_find_edge = -1

    do k=g%ptr(i),g%ptr(i+1)-1
        if (g%node(k)==j) cs_find_edge = k
    enddo

end function cs_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function cs_make_cursor(g,thread) result(cursor)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: k

    cursor%start = 1
    cursor%final = g%capacity
    cursor%current = 0

    k = 1
    do while(g%ptr(k+1)-g%ptr(k)==0)
        k = k+1
    enddo

    cursor%edge = [k, g%node(1)]

end function cs_make_cursor



!--------------------------------------------------------------------------!
subroutine cs_get_edges(g,edges,cursor,num_edges,num_returned)             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(out) :: edges(2,num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i, k

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return; we'll either
    ! return how many edges the user asked for, or, if that amount would
    ! go beyond the final edge that the cursor is allowed to access, all
    ! of the remaining edges
    num_returned = min(num_edges,cursor%final-cursor%current)

    ! Fill the edges array's second row with the edge endpoints
    edges(2,1:num_returned) = &
        & g%node(cursor%current+1:cursor%current+num_returned)

    ! Fill in the edges array's first row with the edge start points
    i = cursor%edge(1)

    do k=1,num_returned
        !! This is going to produce code that cannot be analyzed and
        !! optimized well by the compiler. Would be good to find something
        !! more slick with a predictable access pattern.
        do while(cursor%current >= g%ptr(i+1)-1)
            i = i+1
        enddo

        edges(1,k) = i
        cursor%current = cursor%current+1
    enddo

    cursor%edge = edges(:,num_returned)

end subroutine cs_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k
    logical :: added

    if (.not.g%mutable) then
        print *, 'Attempted to add an edge to an immutable CS graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (.not.g%connected(i,j)) then
        ! Try to see if the new neighbor can be added without
        ! reallocating memory
        added = .false.
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)==0) then
                g%node(k) = j
                added = .true.
                exit
            endif
        enddo

        ! If we successfully added the new edge,
        if (added) then
            ! increment the graph's number of edges
            g%ne = g%ne+1

            ! and update the graphs' maximum degree.
            g%max_degree = max(g%max_degree,k-g%ptr(i)+1)
        else
            ! Otherwise, we need to reallocate the internal structure
            ! of the graph, done in a separate subroutine at the end of
            ! this module.
            call g%add_edge_with_reallocation(i,j)
        endif
    endif

end subroutine cs_add_edge



!--------------------------------------------------------------------------!
subroutine cs_delete_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: jt,k,indx,degree

    if (.not.g%mutable) then
        print *, 'Attempted to delete an edge from an immutable CS graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    ! Find the index in the list of edges of the edge to be removed
    indx = g%find_edge(i,j)

    if (indx/=-1) then
        ! Record the degree of the node from which an edge is to be
        ! removed, and find the last node that i is connected to
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)/=0) then
                degree = k-g%ptr(i)+1
                jt = g%node(k)
            endif
        enddo

        ! If there were more than two edges connected to node i, then
        ! replace the edge (i,j) with the removed edge (i,jt)
        if (degree>1) then
            g%node(indx) = jt
        endif

        ! Remove the last edge (i,jt) connected to i
        g%node( g%ptr(i)+degree-1 ) = 0

        if (degree==g%max_degree) then
            call g%max_degree_update()
        endif

        ! Decrement the number of edges in g
        g%ne = g%ne-1
    endif

end subroutine cs_delete_edge



!--------------------------------------------------------------------------!
subroutine cs_graph_left_permute(g,p,edge_p)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,k,ptr(g%n+1),node(g%capacity)

    ! If the user has asked to see the edge permutation, fill the
    ! first and last constituent arrays
    if (present(edge_p)) then
        allocate(edge_p(3,g%n))
        do i=1,g%n
            edge_p(1,i) = g%ptr(i)
            edge_p(3,i) = g%ptr(i+1)-g%ptr(i)
        enddo
    endif

    ! Find the number of edges for each node under the new ordering
    do i=1,g%n
        ptr(p(i)+1) = g%ptr(i+1)-g%ptr(i)
    enddo

    ! Knowing how many edges each node has in the new ordering, we can
    ! prepare the ptr array
    ptr(1) = 1
    do i=1,g%n
        ptr(i+1) = ptr(i+1)+ptr(i)
    enddo

    ! Shuffle the node array
    do i=1,g%n
        do k=0,g%ptr(i+1)-g%ptr(i)-1
            node( ptr(p(i))+k ) = g%node( g%ptr(i)+k )
        enddo
    enddo

    ! Replace g's versions of ptr and node with the reordered versions
    g%ptr = ptr
    g%node = node

    ! If the user has asked to see the edge permutation, fill the second
    ! constituent array
    if (present(edge_p)) then
        do i=1,g%n
            edge_p(2,i) = g%ptr(p(i))
        enddo
    endif

end subroutine cs_graph_left_permute



!--------------------------------------------------------------------------!
subroutine cs_graph_right_permute(g,p,edge_p)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: j,k

    do k=1,g%capacity
        j = g%node(k)
        if (j/=0) g%node(k) = p(j)
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine cs_graph_right_permute



!--------------------------------------------------------------------------!
subroutine cs_graph_compress(g,edge_p)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, allocatable, intent(inout), optional :: edge_p(:,:)
    ! local variables
    integer :: i,j,k,n
    integer, allocatable :: ptr(:), node(:)

    ! If there are no null nodes stored in g, we don't need to compress
    ! the graph at all
    if (minval(g%node)/=0) then
        if (present(edge_p)) allocate(edge_p(0,0))
    else
        ! Allocate the ptr array
        allocate(ptr(g%n+1),node(g%ne))
        ptr = 0
        ptr(1) = 1
        n = 0

        ! Fill out the ptr array to reflect the number of non-null edges
        ! and copy endpoints of the non-null edges into the node array
        do i=1,g%n
            do k=g%ptr(i),g%ptr(i+1)-1
                j = g%node(k)

                if (j/=0) then
                    ptr(i+1) = ptr(i+1)+1
                    n = n+1
                    node(n) = j
                endif
            enddo

            ptr(i+1) = ptr(i)+ptr(i+1)
        enddo

        if (present(edge_p)) then
            allocate(edge_p(3,g%n))

            ! Fill the edge permutation array
            do i=1,g%n
                edge_p(1,i) = g%ptr(i)
                edge_p(2,i) = ptr(i)
                edge_p(3,i) = ptr(i+1)-ptr(i)
            enddo
        endif

        ! Move the allocation status from the ptr and node arrays to those
        ! owned by the graph
        call move_alloc(from=ptr, to=g%ptr)
        call move_alloc(from=node, to=g%node)

        ! Change the graph's capacity to reflect that there is no
        ! extra storage
        g%capacity = g%ne
    endif

    ! Mark the graph as immutable
    g%mutable = .false.

end subroutine cs_graph_compress



!--------------------------------------------------------------------------!
subroutine cs_graph_decompress(g)                                          !
!--------------------------------------------------------------------------!
    class(cs_graph), intent(inout) :: g

    g%mutable = .true.

end subroutine cs_graph_decompress




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_free(g)                                                      !
!--------------------------------------------------------------------------!
    class(cs_graph), intent(inout) :: g

    deallocate(g%ptr,g%node)
    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0
    g%capacity = 0

    g%mutable = .true.

end subroutine cs_free




!==========================================================================!
!==== Testing, debugging & I/O                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine cs_dump_edges(g,edges)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,k

    do i=1,g%n
        do k=g%ptr(i),g%ptr(i+1)-1
            edges(1,k) = i
            edges(2,k) = g%node(k)
        enddo
    enddo

end subroutine cs_dump_edges




!==========================================================================!
!==== Auxiliary procedures                                             ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine add_edge_with_reallocation(g,i,j)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    ! local variables
    integer :: k
    integer, allocatable :: node(:)

    ! Allocate a temporary array `node` which will store the new list of
    ! destination edges of g
    allocate(node(g%capacity+1))

    ! The index `k` in the new array of the edge to be added is 1 after
    ! the last index for all the edges of row i
    k = g%ptr(i+1)

    ! Copy over the old edges of g into the new array
    node(1:k-1) = g%node(1:k-1)
    node(k+1:g%capacity+1) = g%node(k:g%capacity)

    ! Add in the new edge
    node(k) = j

    ! Transfer the allocation status from the temporary edge array to
    ! g's edge array
    call move_alloc(from=node, to=g%node)

    ! Update all the pointers to the starting index in the array `node`
    ! for the edges of each vertex k
    do k=i+1,g%n+1
        g%ptr(k) = g%ptr(k)+1
    enddo

    ! Increment the number of edges and the capacity of g
    g%ne = g%ne+1
    g%capacity = g%capacity+1

    ! Update the max degree of g if need be
    g%max_degree = max(g%max_degree,g%ptr(i+1)-g%ptr(i))

end subroutine add_edge_with_reallocation



!--------------------------------------------------------------------------!
subroutine sort_node(g)                                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: k,start,finish,num,p(g%max_degree),node(g%max_degree)

    do k=1,g%n
        start  = g%ptr(k)
        finish = g%ptr(k+1)-1
        num = finish-start+1
        node(1:num) = g%node(start:finish)
        p(1:num) = order(node(1:num))
        g%node(start:finish) = node(p(1:num))
    enddo

end subroutine sort_node



!--------------------------------------------------------------------------!
subroutine max_degree_update(g)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(cs_graph), intent(inout) :: g
    ! local variables
    integer :: i,k,degree,max_d

    max_d = 0

    ! loop through all the nodes
    do i=1,g%n
        ! set the degree of i to be 0
        degree = 0

        ! check how many neighbors i has
        do k=g%ptr(i),g%ptr(i+1)-1
            if (g%node(k)/=0) degree = k-g%ptr(i)+1
        enddo

        ! update the max degree of g if node i has a greater degree
        max_d = max(max_d,degree)
    enddo

    g%max_degree = max_d

end subroutine max_degree_update



end module cs_graphs
