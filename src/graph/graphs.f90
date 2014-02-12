!==========================================================================!
!==========================================================================!
module graphs                                                              !
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
type, abstract :: graph                                                    !
!--------------------------------------------------------------------------!
    integer :: n,m,ne,capacity,max_degree
    logical :: mutable
contains
    !--------------
    ! Constructors
    !--------------
    procedure(init_graph_ifc), deferred :: init
    ! Initialize a graph based on the number of left- and right-
    ! vertices, and optionally the maximum number of neighbors for each
    ! vertex.

    procedure(copy_graph_ifc), deferred :: copy
    ! Copy the connectivity structure of another graph, which may well
    ! be of a different type.


    !-----------
    ! Accessors
    !-----------
    procedure(degree_ifc), deferred    :: degree
    ! Return the degree of a given vertex

    procedure(neighbors_ifc), deferred :: neighbors
    ! Return all neighbors of a given vertex

    procedure(connected_ifc), deferred :: connected
    ! Return true if two vertices i, j are connected, false otherwise.

    procedure(find_edge_ifc), deferred :: find_edge
    ! Find the index of the edge between the two vertices (i,j), if it
    ! exists; return -1 if it does not.


    !---------------
    ! Edge iterator
    !---------------
    procedure(make_cursor_ifc), deferred :: make_cursor
    ! Make a cursor which stores some placeholder information needed
    ! for iterating through all of a graph's edges.

    procedure(get_edges_ifc), deferred   :: get_edges
    ! Return a fixed number of edges of the graph and update the graph
    ! edge cursor to reflect our new position with the graph's edges.


    !----------
    ! Mutators
    !----------
    procedure(add_edge_ifc), deferred      :: add_edge
    ! Add in a new edge if it does not already exist.
    ! The behavior of this procedure will change depending on the graph's
    ! mutability state; if the graph has been compressed, it will be
    ! rendered immutable, and attemting to add an edge will instead
    ! yield an error.

    procedure(delete_edge_ifc), deferred   :: delete_edge
    ! Delete an edge if it does exist.
    ! Will result in an error if invoked when the graph is immutable.

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

!    procedure(compress_graph_ifc), deferred :: compress
    ! Compress the graph. This eliminates any extra space which has been
    ! pre-allocated for adding edges. Reduces memory usage and branching
    ! but renders the graph immutable.

!    procedure(add_graph_ifc), deferred      :: add


    !-------------
    ! Destructors
    !-------------
    procedure(free_ifc), deferred :: free
    ! Set all graph attributes to 0 and deallocate any internal data.


    !--------------------------
    ! Testing, debugging & I/O
    !--------------------------
    procedure(dump_edges_ifc), deferred :: dump_edges
    ! Write all of the graph's edges to an array.

    procedure :: write_to_file
    ! Write to a file the number of left- and right-vertices of the
    ! graph, the number of edges, and then all of the edges.

end type graph



!--------------------------------------------------------------------------!
type :: graph_edge_cursor                                                  !
!--------------------------------------------------------------------------!
    integer :: start, final, indx, current, edge(2)
end type graph_edge_cursor



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_graph_ifc(g,n,m,num_neighbor_nodes)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: n
        integer, intent(in), optional :: m, num_neighbor_nodes(:)
    end subroutine init_graph_ifc

    subroutine copy_graph_ifc(g,h)
        import :: graph
        class(graph), intent(inout) :: g
        class(graph), intent(in)    :: h
    end subroutine copy_graph_ifc

    function degree_ifc(g,i) result(d)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i
        integer :: d
    end function degree_ifc

    subroutine neighbors_ifc(g,i,nbrs)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i
        integer, intent(out) :: nbrs(:)
    end subroutine neighbors_ifc

    function connected_ifc(g,i,j)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i,j
        logical :: connected_ifc
    end function connected_ifc

    function find_edge_ifc(g,i,j)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i,j
        integer :: find_edge_ifc
    end function find_edge_ifc

    function make_cursor_ifc(g,thread) result(cursor)
        import :: graph, graph_edge_cursor
        class(graph), intent(in) :: g
        integer, intent(in) :: thread
        type(graph_edge_cursor) :: cursor
    end function make_cursor_ifc

    function get_edges_ifc(g,cursor,num_edges,num_returned) result(edges)
        import :: graph, graph_edge_cursor
        class(graph), intent(in) :: g
        type(graph_edge_cursor), intent(inout) :: cursor
        integer, intent(in) :: num_edges
        integer, intent(out) :: num_returned
        integer :: edges(2,num_edges)
    end function get_edges_ifc

    subroutine add_edge_ifc(g,i,j)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: i,j
    end subroutine add_edge_ifc

    subroutine delete_edge_ifc(g,i,j)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: i,j
    end subroutine delete_edge_ifc

    subroutine permute_graph_ifc(g,p,edge_p)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: p(:)
        integer, allocatable, intent(out), optional :: edge_p(:,:)
    end subroutine permute_graph_ifc

    subroutine compress_graph_ifc(g)
        import :: graph
        class(graph), intent(inout) :: g
    end subroutine compress_graph_ifc

    subroutine add_graph_ifc(g,g1,g2)
        import :: graph
        class(graph), intent(inout) :: g
        class(graph), intent(in)    :: g1, g2
    end subroutine add_graph_ifc

    subroutine free_ifc(g)
        import :: graph
        class(graph), intent(inout) :: g
    end subroutine free_ifc

    subroutine dump_edges_ifc(g,edges)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(out) :: edges(:,:)
    end subroutine dump_edges_ifc
end interface



!--------------------------------------------------------------------------!
type :: graph_pointer                                                      !
!--------------------------------------------------------------------------!
    class(graph), pointer :: g
end type graph_pointer



contains



!--------------------------------------------------------------------------!
subroutine write_to_file(g,filename)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    character(len=*), intent(in) :: filename
    ! local variables
    integer :: n,edges(2,g%ne)

    call g%dump_edges(edges)

    open(unit=10,file=trim(filename))
    write(10,*) g%n, g%ne
    do n=1,g%ne
        write(10,*) edges(:,n)
    enddo
    close(10)

end subroutine write_to_file





end module graphs
