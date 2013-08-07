module wrapper

use iso_c_binding
use graphs
use ll_graphs
use coo_graphs
use cs_graphs



contains



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
subroutine graph_init_c(cgp,n,m) bind(c,name='graph_init')                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: n,m
    ! local variables
    type(graph_pointer), pointer :: gp
    class(graph), pointer :: g

    call c_f_pointer(cgp,gp)
    g => gp%g
    print *, associated(g)
    call g%init(n,m)

end subroutine graph_init_c



!--------------------------------------------------------------------------!
subroutine yea_bitches(cgp) bind(c)                                      !
!--------------------------------------------------------------------------!
    type(c_ptr), intent(in) :: cgp
    type(graph_pointer), pointer :: gp
    class(graph), pointer :: g


    call c_f_pointer(cgp,gp)
    print *, associated(g)
    g => gp%g
    print *, g%n, g%m

end subroutine yea_bitches





end module wrapper
