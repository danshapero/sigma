module ll_graphs

use graph_interface
use types, only: dynamic_array

implicit none


!--------------------------------------------------------------------------!
type, extends(graph) :: ll_graph                                           !
!--------------------------------------------------------------------------!
    type(dynamic_array), allocatable :: lists(:)
contains
    !--------------
    ! Constructors
    procedure :: init_const_degree => ll_init_const_degree
    procedure :: init_variable_degree => ll_init_variable_degree
    procedure :: copy => ll_graph_copy

    !-----------
    ! Accessors
    procedure :: degree => ll_degree
    procedure :: get_neighbors => ll_get_neighbors
    procedure :: connected => ll_connected
    procedure :: find_edge => ll_find_edge

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
    procedure :: compress => ll_graph_compress

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
subroutine ll_init_const_degree(g,n,m,degree)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n
    integer, intent(in), optional :: m, degree
    ! local variables
    integer :: k, ne

    g%n = n
    allocate(g%lists(n))

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! If we know how many neighbors each vertex has, initialize each
    ! dynamic array with the requisite amount of storage.
    if (present(degree)) then
        do k=1,n
            call g%lists(k)%init(capacity=degree, min_capacity=2)
        enddo
        ne = degree*g%n
    else
        do k=1,n
            call g%lists(k)%init(capacity=4, min_capacity=2)
        enddo
        ne = 4*g%n
    endif

    ! The total number of edges and max degree at initialization is zero
    g%ne = 0
    g%max_degree = 0
    g%capacity = ne

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine ll_init_const_degree



!--------------------------------------------------------------------------!
subroutine ll_init_variable_degree(g,n,m,degrees)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: n, degrees(:)
    integer, intent(in), optional :: m
    ! local variables
    integer :: k

    g%n = n
    allocate(g%lists(n))

    if (present(m)) then
        g%m = m
    else
        g%m = n
    endif

    ! Initialize each dynamic array with the storage specified by the
    ! `degrees` argument.
    do k=1,n
        call g%lists(k)%init(capacity=degrees(k), min_capacity=2)
    enddo

    ! The total number of edges and max degree at initialization is zero
    g%ne = 0
    g%max_degree = 0
    g%capacity = sum(degrees)

    ! Mark the graph as mutable
    g%mutable = .true.

end subroutine ll_init_variable_degree



!--------------------------------------------------------------------------!
subroutine ll_graph_copy(g,h,trans)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    class(graph), intent(in)       :: h
    logical, intent(in), optional :: trans
    ! local variables
    integer :: ind(2),order(2),nv(2),k,n,num_returned,num_blocks,edges(2,64)
    type(graph_edge_cursor) :: cursor

    ! Check if we're copying h or h with all directed edges reversed
    nv = [h%n, h%m]
    order = [1, 2]

    if (present(trans)) then
        if (trans) then
            nv = [h%m, h%n]
            order = [2,1]
        endif
    endif

    ! Mark the graph as mutable
    g%mutable = .true.

    ! Initialize g to have the same number of left- and right-nodes as h
    call g%init(nv(1),nv(2))

    ! Get a cursor from h with which to iterate through its edges
    cursor = h%make_cursor(0)

    ! Find the number of chunks into which we're dividing the edges of h
    num_blocks = (cursor%final-cursor%start)/64+1

    ! Iterate through all the chunks
    do n=1,num_blocks
        ! Get a chunk of edges from h
        call h%get_edges(edges,cursor,64,num_returned)

        ! For each edge,
        do k=1,num_returned
            ind = edges(order,k)

            ! If that edge isn't null,
            if (ind(1)/=0 .and. ind(2)/=0) then
                ! Add it to g
                call g%add_edge(ind(1),ind(2))
            endif
        enddo
    enddo

end subroutine ll_graph_copy




!==========================================================================!
!==== Accessors                                                        ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ll_degree(g,i) result(d)                                          !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i
    integer :: d

    d = g%lists(i)%length

end function ll_degree



!--------------------------------------------------------------------------!
subroutine ll_get_neighbors(g,neighbors,i)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: neighbors(:)
    integer, intent(in) :: i
    ! local variables
    integer :: k

    neighbors = 0
    do k=1,g%lists(i)%length
        neighbors(k) = g%lists(i)%get_entry(k)
    enddo

end subroutine ll_get_neighbors



!--------------------------------------------------------------------------!
function ll_connected(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i,j
    logical :: ll_connected
    ! local variables
    integer :: k

    ll_connected = .false.
    do k=1,g%lists(i)%length
        if (g%lists(i)%get_entry(k)==j) then
            ll_connected = .true.
            exit
        endif
    enddo

end function ll_connected



!--------------------------------------------------------------------------!
function ll_find_edge(g,i,j)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: i,j
    integer :: ll_find_edge
    ! local variables
    integer :: k,total

    ll_find_edge = -1

    total = 0
    do k=1,i-1
        total = total+g%lists(k)%length
    enddo

    do k=1,g%lists(i)%length
        if (g%lists(i)%get_entry(k)==j) then
            ll_find_edge = total+k
            exit
        endif
    enddo

end function ll_find_edge




!==========================================================================!
!==== Edge iterator                                                    ====!
!==========================================================================!

!--------------------------------------------------------------------------!
function ll_make_cursor(g,thread) result(cursor)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(in) :: thread
    type(graph_edge_cursor) :: cursor
    ! local variables
    integer :: k

    cursor%start = 1
    cursor%final = g%ne
    cursor%current = 0

    k = 1
    do while (g%lists(k)%length==0)
        k = k+1
    enddo

    cursor%edge = [k, g%lists(k)%get_entry(1)]
    cursor%indx = 0

end function ll_make_cursor



!--------------------------------------------------------------------------!
subroutine ll_get_edges(g,edges,cursor,num_edges,num_returned)             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: edges(2,num_edges)
    type(graph_edge_cursor), intent(inout) :: cursor
    integer, intent(in) :: num_edges
    integer, intent(out) :: num_returned
    ! local variables
    integer :: i,k,num_added,num_from_this_row

    ! Set up the returned edges to be 0
    edges = 0

    ! Count how many edges we're actually going to return
    num_returned = min(num_edges, cursor%final-cursor%current)

    ! Find the last row we left off at
    i = cursor%edge(1)

    num_added = 0
    do while(num_added<num_returned)
        ! Either we are returning all the edges for the current node i,
        ! or the current request doesn't call for so many edges and we are
        ! returning fewer
        num_from_this_row = min(g%lists(i)%length-cursor%indx, &
                                & num_returned-num_added)

        do k=1,num_from_this_row
            edges(1,num_added+k) = i
            edges(2,num_added+k) = g%lists(i)%get_entry(cursor%indx+k)
        enddo

        !! Check that this is right
        !cursor%indx = mod(cursor%indx+num_from_this_row,g%lists(i)%length)

        ! If we returned all nodes from this row, increment the row
        if (num_from_this_row == g%lists(i)%length-cursor%indx) then
            i = i+1
            cursor%indx = 0
        else
            cursor%indx = cursor%indx+num_from_this_row
        endif

        ! Increase the number of edges added
        num_added = num_added+num_from_this_row
    enddo

    cursor%current = cursor%current+num_returned
    cursor%edge(1) = i

end subroutine ll_get_edges




!==========================================================================!
!==== Mutators                                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_add_edge(g,i,j)                                              !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    integer :: cap

    if (.not.g%mutable) then
        print *, 'Attempted to add an edge to an immutable LL graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (.not.g%connected(i,j)) then
        ! Record the old capacity of the list of i's neighbors
        cap = g%lists(i)%capacity

        ! Push the new neighbor node j onto the list of i's neighbors
        call g%lists(i)%push(j)

        ! Change the graph's max degree if need be
        g%max_degree = max(g%max_degree,g%lists(i)%length)

        ! Increment the number of edges and graph capacity
        g%ne = g%ne+1
        g%capacity = g%capacity+g%lists(i)%capacity-cap
    endif

end subroutine ll_add_edge



!--------------------------------------------------------------------------!
subroutine ll_delete_edge(g,i,j)                                           !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: i,j
    integer :: k,jt,degree,cap

    if (.not.g%mutable) then
        print *, 'Attempted to delete an edge from an immutable LL graph'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (g%connected(i,j)) then
        ! Record the degree of vertex i and capacity
        degree = g%lists(i)%length
        cap = g%lists(i)%capacity

        ! Pop from the list of i's neighbors
        jt = g%lists(i)%pop()

        ! If the vertex jt popped from i's neighbors is not vertex j,
        if (jt/=j) then
            ! find where vertex j was stored and put jt there.
            do k=1,g%lists(i)%length
                if (g%lists(i)%get_entry(k)==j) then
                    call g%lists(i)%set_entry(k,jt)
                endif
            enddo
        endif

        ! If the degree of vertex i was the max degree of the graph,
        ! check that the max degree of the graph hasn't decreased.
        !!Make this a guaranteed O(1) operation somehow
        if (degree==g%max_degree) then
            g%max_degree = 0
            do k=1,g%n
                g%max_degree = max(g%max_degree,g%lists(k)%length)
            enddo
        endif

        ! Decrement the number of edges and capacity of the graph
        g%ne = g%ne-1
        g%capacity = g%capacity+g%lists(i)%capacity-cap
    endif

end subroutine ll_delete_edge



!--------------------------------------------------------------------------!
subroutine ll_graph_left_permute(g,p,edge_p)                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,j,k
    type(dynamic_array) :: lists(g%n)

    ! If the user wants to see the edge permutation, set up the first
    ! constituent array before performing the permutation
    if (present(edge_p)) then
        allocate(edge_p(3,g%n))
        edge_p(1,1) = 1
        do i=1,g%n-1
            edge_p(1,i+1) = edge_p(1,i)+g%lists(i)%length
        enddo
    endif

    ! Copy all of the structure of g to a new list of lists and empty g
    do i=1,g%n
        call lists(i)%init(capacity=g%lists(i)%capacity,min_capacity=2)
        do k=1,g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            call lists(i)%push(j)
        enddo
        call g%lists(i)%free()
    enddo

    ! Refill g with the appropriately permuted structure
    do i=1,g%n
        call g%lists(p(i))%init(capacity=lists(i)%capacity,min_capacity=2)
        do k=1,lists(i)%length
            j = lists(i)%get_entry(k)
            call g%lists(p(i))%push(j)
        enddo
        call lists(i)%free()
    enddo

    ! If the user wants to see the edge permutation, finish the remaining
    ! constituent arrays
    if (present(edge_p)) then
        edge_p(3,1) = 1
        do i=1,g%n-1
            edge_p(3,i+1) = edge_p(3,i)+g%lists(i)%length
        enddo

        do i=1,g%n
            edge_p(2,i) = edge_p(3,p(i))
        enddo

        do i=1,g%n
            edge_p(3,i) = g%lists(p(i))%length
        enddo
    endif

end subroutine ll_graph_left_permute



!--------------------------------------------------------------------------!
subroutine ll_graph_right_permute(g,p,edge_p)                              !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(inout) :: g
    integer, intent(in) :: p(:)
    integer, allocatable, intent(out), optional :: edge_p(:,:)
    ! local variables
    integer :: i,j,k

    do i=1,g%n
        do k=1,g%lists(i)%length
            j = g%lists(i)%get_entry(k)
            if (j/=0) call g%lists(i)%set_entry(k,p(j))
        enddo
    enddo

    if (present(edge_p)) allocate(edge_p(0,0))

end subroutine ll_graph_right_permute



!--------------------------------------------------------------------------!
subroutine ll_graph_compress(g,edge_p)                                     !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer, allocatable, intent(inout), optional :: edge_p(:,:)

    if (present(edge_p)) allocate(edge_p(0,0))

    ! LL graphs cannot have their storage compressed.

    ! Mark the graph as immutable
    g%mutable = .false.

end subroutine ll_graph_compress




!==========================================================================!
!==== Destructors                                                      ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_destroy(g)                                                   !
!--------------------------------------------------------------------------!
    class(ll_graph), intent(inout) :: g
    integer :: i

    do i=1,g%n
        call g%lists(i)%free()
    enddo

    deallocate(g%lists)

    g%n = 0
    g%m = 0
    g%ne = 0
    g%max_degree = 0
    g%capacity = 0

    g%mutable = .true.

end subroutine ll_destroy




!==========================================================================!
!==== Testing, debugging & I/O                                         ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine ll_dump_edges(g,edges)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(ll_graph), intent(in) :: g
    integer, intent(out) :: edges(:,:)
    ! local variables
    integer :: i,k,next

    next = 0
    do i=1,g%n
        do k=1,g%lists(i)%length
            next = next+1
            edges(1,next) = i
            edges(2,next) = g%lists(i)%get_entry(k)
        enddo
    enddo

end subroutine ll_dump_edges






end module ll_graphs
