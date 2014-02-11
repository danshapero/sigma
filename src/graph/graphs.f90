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
    ! core graph procedures
    procedure(init_graph_ifc), deferred     :: init
    procedure(copy_graph_ifc), deferred     :: copy
    procedure(degree_ifc), deferred         :: degree
    procedure(neighbors_ifc), deferred      :: neighbors
    procedure(connected_ifc), deferred      :: connected
    procedure(find_edge_ifc), deferred      :: find_edge
    procedure(make_cursor_ifc), deferred    :: make_cursor
    procedure(get_edges_ifc), deferred      :: get_edges
    procedure(add_edge_ifc), deferred       :: add_edge
    procedure(delete_edge_ifc), deferred    :: delete_edge
    procedure(permute_graph_ifc), deferred  :: left_permute
    procedure(permute_graph_ifc), deferred  :: right_permute
!    procedure(compress_graph_ifc), deferred :: compress
!    procedure(add_graph_ifc), deferred      :: add
    procedure(free_ifc), deferred           :: free
    procedure(dump_edges_ifc), deferred     :: dump_edges
    ! procedures for testing, debugging and i/o
    procedure :: write_to_file
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
