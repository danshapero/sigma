module sigma

! Use all the graph modules
use graph_interface
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs
use graphs

! Use all the linear operator modules
use linear_operator_interface
use linear_operator_sums
use linear_operator_products
use linear_operators

! Use all the matrix modules
use sparse_matrices

! Use the solver and preconditioner modules
use cg_solvers
use bicgstab_solvers
use jacobi_solvers

! Use the C wrapper module
use wrapper

! Use other auxiliary modules
use eigensolver
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
