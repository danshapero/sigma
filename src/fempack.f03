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



contains


!--------------------------------------------------------------------------!
subroutine new_graph(g,graph_format,n,m,edges)                             !
!--------------------------------------------------------------------------!
    class(graph), pointer, intent(inout) :: g
    character(len=*), intent(in) :: graph_format
    integer, intent(in) :: n
    integer, intent(in), optional :: m, edges(:,:)

    select case(trim(graph_format))
        case('cs')
            allocate(cs_graph::g)
        case('coo')
            allocate(coo_graph::g)
        case('ll')
            allocate(ll_graph::g)
    end select

    call g%init(n,m,edges)
    
end subroutine new_graph



!--------------------------------------------------------------------------!
subroutine new_sparse_matrix(A,matrix_format,nrow,ncol,g)                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), pointer, intent(inout) :: A
    character(len=*), intent(in) :: matrix_format
    integer, intent(in) :: nrow, ncol
    class(graph), pointer, intent(in), optional :: g
    ! local variables
    class(graph), pointer :: ag

    select case(trim(matrix_format))
        case('csr')
            allocate(cs_matrix::A)
            A%orientation = 'row'
        case('csc')
            allocate(cs_matrix::A)
            A%orientation = 'col'
        case('coo')
            allocate(coo_matrix::A)
            A%orientation = 'row'
    end select

    A%nrow = nrow
    A%ncol = ncol

    if (present(g)) then
        call A%assemble(g)
    else
        select type(A)
            type is(cs_matrix)
                allocate(cs_graph::ag)
            type is(coo_matrix)
                allocate(coo_graph::ag)
        end select

        select case(A%orientation)
            case('row')
                call ag%init(nrow,ncol)
            case('col')
                call ag%init(ncol,nrow)
        end select

        call A%assemble(ag)

        ! nullify(ag)   <- should we do this?
    endif

end subroutine new_sparse_matrix




! Write this one too:
!function new_block_matrix



end module fempack

