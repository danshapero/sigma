module wrapper

use iso_c_binding
use graphs
use ll_graphs
use coo_graphs
use cs_graphs



contains



!==========================================================================!
!==========================================================================!
!==== Routines for C <-> F2K pointer translation                       ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine get_graph(cgp,storage_format) bind(c)                           !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(inout) :: cgp
    integer(c_int), intent(in), value :: storage_format
    ! local variables
    type(graph_pointer), pointer :: gp
    !class(graph), pointer :: g

    allocate(gp)
    cgp = c_loc(gp)

    select case(storage_format)
        case(0)
            allocate(ll_graph  ::gp%g)
        case(1)
            allocate(coo_graph ::gp%g)
        case(2)
            allocate(cs_graph  ::gp%g)
    end select

end subroutine get_graph



!--------------------------------------------------------------------------!
function cgraph_to_fgraph(cgp) result(g)                                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    class(graph), pointer :: g
    ! local variables
    type(graph_pointer), pointer :: gp

    call c_f_pointer(cgp,gp)
    g => gp%g

end function cgraph_to_fgraph








!==========================================================================!
!==========================================================================!
!==== C-compatible wrappers to Fortran graph operations                ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine graph_init_c(cgp,n,m) bind(c,name='graph_init')                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: n,m
    ! local variables
    class(graph), pointer :: g

    g => cgraph_to_fgraph(cgp)
    call g%init(n,m)

end subroutine graph_init_c



!--------------------------------------------------------------------------!
subroutine connected_c(cgp,i,j,con) bind(c,name='connected')               !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value  :: i,j
    integer(c_int), intent(out) :: con
    ! local variables
    class(graph), pointer :: g

    g => cgraph_to_fgraph(cgp)
    if (g%connected(i,j)) then
        con = 1
    else
        con = 0
    endif

end subroutine connected_c



!--------------------------------------------------------------------------!
subroutine add_edge_c(cgp,i,j) bind(c,name='add_edge')                     !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: i,j
    ! local variables
    class(graph), pointer :: g

    g => cgraph_to_fgraph(cgp)
    call g%add_edge(i,j)

end subroutine add_edge_c



!--------------------------------------------------------------------------!
subroutine delete_edge_c(cgp,i,j) bind(c,name='delete_edge')               !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: i,j
    ! local variables
    class(graph), pointer :: g

    g => cgraph_to_fgraph(cgp)
    call g%delete_edge(i,j)

end subroutine delete_edge_c







end module wrapper
