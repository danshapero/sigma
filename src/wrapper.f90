module wrapper

use iso_c_binding

use graphs
use ll_graphs
use coo_graphs
use cs_graphs
use ellpack_graphs

use sparse_matrices
use coo_matrices
use cs_matrices
use ellpack_matrices

use conversions


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

    allocate(gp)
    cgp = c_loc(gp)

    select case(storage_format)
        case(0)
            allocate(ll_graph::gp%g)
        case(1)
            allocate(coo_graph::gp%g)
        case(2)
            allocate(cs_graph::gp%g)
        case(3)
            allocate(ellpack_graph::gp%g)
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



!--------------------------------------------------------------------------!
subroutine get_matrix(cmp,storage_format) bind(c)                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(inout) :: cmp
    integer(c_int), intent(in), value :: storage_format
    ! local variables
    type(sparse_matrix_pointer), pointer :: mp

    allocate(mp)
    cmp = c_loc(mp)

    select case(storage_format)
        case(0)
            allocate(coo_matrix::mp%A)
        case(1)
            allocate(cs_matrix::mp%A)
        case(2)
            allocate(ellpack_matrix::mp%A)
    end select

end subroutine get_matrix



!--------------------------------------------------------------------------!
function cmatrix_to_fmatrix(cmp) result(A)                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cmp
    class(sparse_matrix), pointer :: A
    ! local variables
    type(sparse_matrix_pointer), pointer :: mp

    call c_f_pointer(cmp,mp)
    A => mp%A

end function cmatrix_to_fmatrix








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
subroutine graph_neighbors_c(cgp,i,nbrs,d) bind(c,name='graph_neighbors')  !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: i,d
    integer(c_int), intent(out) :: nbrs(d)
    ! local variables
    class(graph), pointer :: g

    g => cgraph_to_fgraph(cgp)
    call g%neighbors(i+1,nbrs)
    nbrs = nbrs-1

end subroutine graph_neighbors_c



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
    if (g%connected(i+1,j+1)) then
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
    call g%add_edge(i+1,j+1)

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
    call g%delete_edge(i+1,j+1)

end subroutine delete_edge_c



!--------------------------------------------------------------------------!
subroutine convert_c(cgp,storage_format) bind(c,name='convert')            !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cgp
    integer(c_int), intent(in), value :: storage_format
    ! local variables
    type(graph_pointer), pointer :: gp
    character(len=3) :: str_fmt

    select case(storage_format)
        case(0)
            str_fmt = 'll '
        case(1)
            str_fmt = 'coo'
        case(2)
            str_fmt = 'cs '
    end select

    call c_f_pointer(cgp,gp)
    call convert(gp%g,str_fmt)

end subroutine convert_c






!==========================================================================!
!==========================================================================!
!==== C-compatible wrappers to Fortran matrix routines                 ====!
!==========================================================================!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine sparse_matrix_init_c(cmp,nrow,ncol,orientation) &               !
    & bind(c,name='sparse_matrix_init')                                    !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cmp
    integer(c_int), intent(in), value :: nrow,ncol
    integer(c_int), intent(in), value :: orientation
    ! local variables
    class(sparse_matrix), pointer :: A
!    class(graph), pointer :: g

    A => cmatrix_to_fmatrix(cmp)
!    g => cgraph_to_fgraph(cgp)
    select case(orientation)
        case(1)
            call A%init(nrow,ncol,'row')
        case(-1)
            call A%init(nrow,ncol,'col')
    end select

end subroutine sparse_matrix_init_c



!--------------------------------------------------------------------------!
subroutine sparse_matrix_set_value_c(cmp,i,j,val) &                        !
    & bind(c,name='sparse_matrix_set_value')                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    type(c_ptr), intent(in) :: cmp
    integer(c_int), intent(in), value :: i,j
    real(c_double), intent(in), value :: val
    ! local variables
    class(sparse_matrix), pointer :: A

    A => cmatrix_to_fmatrix(cmp)

    call A%set_value(i,j,val)

end subroutine sparse_matrix_set_value_c




end module wrapper
