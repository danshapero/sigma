!==========================================================================!
!==========================================================================!
module graph_interfaces                                                    !
!==========================================================================!
!==========================================================================!
!==== This module contains the definition of the abstract graph data   ====!
!==== type, which is used as one of the key underlying data structures ====!
!==== of sparse matrices. The graph data type is abstract, so this     ====!
!==== module describes only the interface for graph objects; the       ====!
!==== implementations of the graph interface are contained in all the  ====!
!==== other files in this directory, e.g. cs_graphs.f90, etc.          ====!
!==========================================================================!
!==========================================================================!

implicit none




!--------------------------------------------------------------------------!
type, abstract :: graph_interface                                          !
!--------------------------------------------------------------------------!
    integer :: n, m, reference_count = 0
contains
    !--------------
    ! Constructors
    !--------------
    procedure(init_graph_ifc), deferred :: init
    ! Initialize an empty graph

    procedure(copy_graph_ifc), deferred :: copy
    ! Copy the connectivity structure of another graph, which may well
    ! be of a different type.


    !-----------
    ! Accessors
    !-----------
    procedure(get_num_edges_ifc), deferred :: get_num_edges
    ! Return the total number of edges of the graph

    procedure(get_degree_ifc), deferred :: get_degree
    ! Return the degree of a given vertex

    procedure(get_num_edges_ifc), deferred :: get_max_degree
    ! Return the maximum degree among all vertices of the graph

    procedure(get_neighbors_ifc), deferred :: get_neighbors
    ! Return all neighbors of a given vertex

    procedure(connected_ifc), deferred :: connected
    ! Return true if two vertices i, j are connected, false otherwise.

    procedure(find_edge_ifc), deferred :: find_edge
    ! Find the index of the edge between the two vertices (i,j), if it
    ! exists; return -1 if it does not.

    procedure, nopass :: is_get_neighbors_fast => get_neighbors_is_not_fast
    ! Returns true if the graph is in a storage format for which getting
    ! all the neighbors of a vertex can be done in O(max_degree) time;
    ! used for optimizing sparse matrix-matrix multiplication


    !---------------
    ! Edge iterator
    !---------------
    procedure(make_cursor_ifc), deferred :: make_cursor
    ! Make a cursor which stores some placeholder information needed
    ! for iterating through all of a graph's edges.

    procedure(get_edges_ifc), deferred :: get_edges
    ! Return a fixed number of edges of the graph and update the graph
    ! edge cursor to reflect our new position with the graph's edges.


    !----------
    ! Mutators
    !----------
    procedure(change_edge_ifc), deferred :: add_edge
    ! Add in a new edge if it does not already exist.

    procedure(change_edge_ifc), deferred :: delete_edge
    ! Delete an edge if it does exist.

    procedure(permute_graph_ifc), deferred :: left_permute
    ! Apply a permutation to the graph's left-vertices.
    ! Optionally return an array which describes, in compact form, the
    ! resulting permutation to the graph's edges. This optional argument
    ! is necessary when permuting matrices, which need to know how to
    ! rearrange their non-zero entries after changing the underlying
    ! structure.

    procedure(permute_graph_ifc), deferred :: right_permute
    ! Apply a permutation to the graph's right-vertices.
    ! Optionally return compact array describing edge permutation.


    !------------------------------------
    ! Destructors and reference counting
    !------------------------------------
    procedure :: add_reference => graph_add_reference
    ! Whenever another object, such as a matrix, points to a graph, we
    ! need to increment the graph's reference counter.
    !TODO If the reference
    ! counter goes above 2, the graph needs to be made immutable.

    procedure :: remove_reference => graph_remove_reference
    ! If another object that was pointing to a graph is destroyed, the
    ! graph's reference counter decreases.

    procedure(destroy_graph_ifc), deferred :: destroy
    ! Set all graph attributes to 0 and deallocate any internal data.


    !--------------------------
    ! Testing, debugging & I/O
    !--------------------------
    procedure(dump_edges_ifc), deferred :: dump_edges
    ! Write all of the graph's edges to an array.

    procedure :: to_dense_graph
    ! Convert the graph to a dense array.

    procedure :: write_to_file => write_graph_to_file
    ! Write to a file the number of left- and right-vertices of the
    ! graph, the number of edges, and then all of the edges.

end type graph_interface



!--------------------------------------------------------------------------!
type :: graph_edge_cursor                                                  !
!--------------------------------------------------------------------------!
    integer :: first, last, current, idx, edge(2)
    type(graph_edge_cursor), pointer :: sub_cursor(:,:)
contains
    procedure :: done => graph_edge_cursor_is_done
end type graph_edge_cursor



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_graph_ifc(g, n, m)
        import :: graph_interface
        class(graph_interface), intent(inout) :: g
        integer, intent(in) :: n
        integer, intent(in), optional :: m
    end subroutine init_graph_ifc

    subroutine copy_graph_ifc(g, h, trans)
        import :: graph_interface
        class(graph_interface), intent(inout) :: g
        class(graph_interface), intent(in)    :: h
        logical, intent(in), optional :: trans
    end subroutine copy_graph_ifc

    function get_num_edges_ifc(g) result(k)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer :: k
    end function get_num_edges_ifc

    function get_degree_ifc(g, i) result(d)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer, intent(in) :: i
        integer :: d
    end function get_degree_ifc

    subroutine get_neighbors_ifc(g, neighbors, i)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer, intent(in) :: i
        integer, intent(out) :: neighbors(:)
    end subroutine get_neighbors_ifc

    function connected_ifc(g, i, j)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer, intent(in) :: i, j
        logical :: connected_ifc
    end function connected_ifc

    function find_edge_ifc(g, i, j)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer, intent(in) :: i, j
        integer :: find_edge_ifc
    end function find_edge_ifc

    function make_cursor_ifc(g) result(cursor)
        import :: graph_interface, graph_edge_cursor
        class(graph_interface), intent(in) :: g
        type(graph_edge_cursor) :: cursor
    end function make_cursor_ifc

    subroutine get_edges_ifc(g, edges, cursor, num_edges, num_returned)
        import :: graph_interface, graph_edge_cursor
        class(graph_interface), intent(in) :: g
        integer, intent(in) :: num_edges
        integer, intent(out) :: edges(2, num_edges)
        type(graph_edge_cursor), intent(inout) :: cursor
        integer, intent(out) :: num_returned
    end subroutine get_edges_ifc

    subroutine change_edge_ifc(g, i, j)
        import :: graph_interface
        class(graph_interface), intent(inout) :: g
        integer, intent(in) :: i, j
    end subroutine change_edge_ifc

    subroutine permute_graph_ifc(g, p, edge_p)
        import :: graph_interface
        class(graph_interface), intent(inout) :: g
        integer, intent(in) :: p(:)
        integer, allocatable, intent(out), optional :: edge_p(:,:)
    end subroutine permute_graph_ifc

    subroutine destroy_graph_ifc(g)
        import :: graph_interface
        class(graph_interface), intent(inout) :: g
    end subroutine destroy_graph_ifc

    subroutine dump_edges_ifc(g, edges)
        import :: graph_interface
        class(graph_interface), intent(in) :: g
        integer, intent(out) :: edges(:,:)
    end subroutine dump_edges_ifc
end interface



! Parameter specifying the batch size 
integer, parameter :: batch_size = 64



contains



!--------------------------------------------------------------------------!
function get_neighbors_is_not_fast() result(fast)                          !
!--------------------------------------------------------------------------!
    logical :: fast

    fast = .false.

end function get_neighbors_is_not_fast



!--------------------------------------------------------------------------!
function get_neighbors_is_fast() result(fast)                              !
!--------------------------------------------------------------------------!
    logical :: fast

    fast = .true.

end function get_neighbors_is_fast



!--------------------------------------------------------------------------!
subroutine graph_add_reference(g)                                          !
!--------------------------------------------------------------------------!
    class(graph_interface), intent(inout) :: g

    g%reference_count = g%reference_count + 1
    !TODO reimplement: if (g%reference_count > 1) g%mutable = .false.

end subroutine graph_add_reference



!--------------------------------------------------------------------------!
subroutine graph_remove_reference(g)                                       !
!--------------------------------------------------------------------------!
    class(graph_interface), intent(inout) :: g

    g%reference_count = g%reference_count - 1

end subroutine graph_remove_reference



!--------------------------------------------------------------------------!
subroutine to_dense_graph(g, A, trans)                                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(in) :: g
    integer, intent(out) :: A(:,:)
    logical, intent(in), optional :: trans
    ! local variables
    integer :: i, j, k, ord(2)
    integer :: num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor

    ord = [1,2]
    if (present(trans)) then
        if (trans) ord = [2,1]
    endif

    ! Set the dense array to 0
    A = 0

    cursor = g%make_cursor()

    ! Iterate through all the edges (i, j) of g
    do while(.not. cursor%done())
        call g%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(ord(1), k)
            j = edges(ord(2), k)

            A(i, j) = 1
        enddo
    enddo

end subroutine to_dense_graph



!--------------------------------------------------------------------------!
subroutine write_graph_to_file(g, filename)                                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph_interface), intent(in) :: g
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: i, j, k
    type(graph_edge_cursor) :: cursor
    integer :: num_returned, edges(2, batch_size)

    !TODO make sure this unit number is safe
    open(unit=10,file=trim(filename))
    write(10,*) g%n, g%m, g%get_num_edges()

    cursor = g%make_cursor()

    do while(.not. cursor%done())
        call g%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            write(10,*) i, j
        enddo
    enddo

    close(10)

end subroutine write_graph_to_file



!--------------------------------------------------------------------------!
pure function graph_edge_cursor_is_done(cursor) result(done)               !
!--------------------------------------------------------------------------!
    class(graph_edge_cursor), intent(in) :: cursor
    logical :: done

    done = cursor%current == cursor%last

end function graph_edge_cursor_is_done





end module graph_interfaces

