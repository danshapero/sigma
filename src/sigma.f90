module sigma

! Use all the graph modules
use graphs
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs

! Use all the matrix modules
use sparse_matrices

! Use the solver and preconditioner modules
use iterative_solvers
use cg_solvers
use bicgstab_solvers
use jacobi_preconditioners

! Use the C wrapper module
use wrapper

! Use other auxiliary modules
use permutations
use conversions
use meshes
use vectors


implicit none



contains


!--------------------------------------------------------------------------!
subroutine new_graph(g,graph_format,n,m,num_nbrs)                          !
!--------------------------------------------------------------------------!
    class(graph), pointer, intent(inout) :: g
    character(len=*), intent(in) :: graph_format
    integer, intent(in) :: n
    integer, intent(in), optional :: m, num_nbrs(:)

    select case(trim(graph_format))
        case('cs')
            allocate(cs_graph::g)
        case('coo')
            allocate(coo_graph::g)
        case('ll')
            allocate(ll_graph::g)
        case('ellpack')
            allocate(ellpack_graph::g)
    end select

    call g%init(n,m=m,degrees=num_nbrs)
    
end subroutine new_graph




end module sigma
