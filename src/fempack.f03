module fempack

! Use all the graph modules
use graphs
use ll_graphs
use coo_graphs
use cs_graphs

! Use all the matrix modules
use sparse_matrices
use cs_matrices
use coo_matrices

! Use all the block matrix modules
use block_sparse_matrices
use bcsr_matrices

! Use the solver and preconditioner modules
use iterative_solvers
use cg_solvers
use jacobi_preconditioners

! Use the C wrapper module
use wrapper

! Use other auxiliary modules
use conversions
use meshes


implicit none



!--------------------------------------------------------------------------!
function new_graph(graph_format,n,m,edges) return(g)                       !
!--------------------------------------------------------------------------!
    character(len=*), intent(in) :: graph_format
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)
    class(graph), allocatable :: g

    select case(trim(graph_format))
        case('cs')
            allocate(cs_graph::g)
        case('coo')
            allocate(coo_graph::g)
        case('ll')
            allocate(ll_graph::g)
    end select

    call g%init(n,m,edges)
    
end function new_graph



!--------------------------------------------------------------------------!
function new_sparse_matrix(matrix_format,nrow,ncol,g) return(A)            !
!--------------------------------------------------------------------------!
    ! input/output variables
    character(len=*), intent(in) :: matrix_format
    integer, intent(in) :: nrow, ncol
    class(graph), pointer, intent(in), optional :: g
    class(sparse_matrix), allocatable :: A

    select case(trim(matrix_format))
        case('csr')
            allocate(cs_matrix::A)
            A%orientation = 'r'
        case('csc')
            allocate(cs_matrix::A)
            A%orientation = 'c'
        case('coo')
            allocate(coo_matrix::A
    end select

end function new_sparse_matrix




! Write this one too:
!function new_block_matrix



end module fempack
