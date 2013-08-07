module graphs

implicit none



!--------------------------------------------------------------------------!
type, abstract :: graph                                                    !
!--------------------------------------------------------------------------!
    integer :: n,m,ne,max_degree
    contains
    procedure(init_graph_ifc), deferred     :: init
    procedure(neighbors_ifc), deferred      :: neighbors
    procedure(connected_ifc), deferred      :: connected
    procedure(find_edge), deferred          :: find_edge
    procedure(add_edge_ifc), deferred       :: add_edge
    procedure(delete_edge_ifc), deferred    :: delete_edge
    procedure(free_ifc), deferred           :: free
    procedure(dump_edges_ifc), deferred     :: dump_edges
end type graph


!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine init_graph_ifc(g,n,m,edges)
        import :: graph
        class(graph), intent(inout) :: g
        integer, intent(in) :: n
        integer, intent(in), optional :: m, edges(:,:)
    end subroutine init_graph_ifc

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

    function find_edge(g,i,j)
        import :: graph
        class(graph), intent(in) :: g
        integer, intent(in) :: i,j
        integer :: find_edge
    end function find_edge

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





end module graphs
