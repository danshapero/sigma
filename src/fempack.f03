module fempack

use graphs
use ll_graphs
use coo_graphs
use cs_graphs

use sparse_matrices
use csr_matrices

implicit none


contains


!--------------------------------------------------------------------------!
subroutine convert(g,storage_format)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), pointer, intent(inout) :: g
    character(len=3), intent(in), optional :: storage_format
    ! local variables
    class(graph), pointer :: gc
    integer :: edges(2,g%ne)
    character(len=3) :: str_fmt

    if (present(storage_format)) then
        str_fmt = storage_format
    else
        str_fmt = 'cs '
    endif

    gc => g
    nullify(g)

    call gc%dump_edges(edges)
    select case(trim(str_fmt))
        case('ll')
            allocate(ll_graph::g)
        case('cs')
            allocate(cs_graph::g)
        case('coo')
            allocate(coo_graph::g)
    end select

    call g%init(gc%n,gc%m,edges)

    call gc%free()
    deallocate(gc)

end subroutine convert



end module fempack
