!==========================================================================!
!==========================================================================!
module default_sparse_matrix_kernels                                       !
!==========================================================================!
!==========================================================================!
!====     This module contains implementations of sparse matrix        ====!
!==== operations which work using only methods provided by the graph   ====!
!==== interface. These computational kernels are used by the default   ====!
!==== sparse matrix implementation, but also by some specific sparse   ====!
!==== matrix classes for which there may be either no performance gain ====!
!==== to be had for a specific kernel, or for which I am too lazy or   ====!
!==== stupid to write that kernel should it exist.                     ====!
!==========================================================================!
!==========================================================================!


use types, only: dp
use graph_interface

implicit none



!--------------------------------------------------------------------------!
abstract interface                                                         !
!--------------------------------------------------------------------------!
    subroutine get_slice_kernel(g, val, nodes, slice, k)
        import :: graph, dp
        class(graph), intent(in) :: g
        real(dp), intent(in) :: val(:)
        integer, intent(out) :: nodes(:)
        real(dp), intent(out) :: slice(:)
        integer, intent(in) :: k
    end subroutine get_slice_kernel

    subroutine permute_kernel(g, val, p)
        import :: graph, dp
        class(graph), intent(inout) :: g
        real(dp), intent(inout) :: val(:)
        integer, intent(in) :: p(:)
    end subroutine permute_kernel

    subroutine matvec_kernel(g, val, x, y)
        import :: graph, dp
        class(graph), intent(in) :: g
        real(dp), intent(in) :: val(:)
        real(dp), intent(in)    :: x(:)
        real(dp), intent(inout) :: y(:)
    end subroutine matvec_kernel
end interface




contains



!==========================================================================!
!==== Accessor kernels                                                 ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine get_slice_contiguous(g, val, nodes, slice, k)                   !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    real(dp), intent(in) :: val(:)
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: l, d, ind

    ! Set the return values to 0
    vals = 0.0_dp

    ! Get the degree of node k
    d = g%degree(k)

    ! Get all the neighbors of k and put them into the array `nodes`
    call g%get_neighbors(nodes, k)

    ! For each neighbor l of node k,
    do ind = 1, d
        l = nodes(ind)

        ! put the matrix entry A(k,l) into the array `slice`
        slice(ind) = val( g%find_edge(k, l) )
    enddo

end subroutine get_slice_contiguous



!--------------------------------------------------------------------------!
subroutine get_slice_discontiguous(g, val, nodes, slice, k)                !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    real(dp), intent(in) :: val(:)
    integer, intent(out) :: nodes(:)
    real(dp), intent(out) :: slice(:)
    integer, intent(in) :: k
    ! local variables
    integer :: l, ind, next

    ! Set the index in `nodes` and `slice` of the next entry to 0
    next = 0

    ! Set the return values to 0
    slice = 0.0_dp

    ! For each node l,
    do l = 1, g%n
        ! Check whether (l, k) is an edge of g
        ind = g%find_edge(l, k)

        ! If it is,
        if (ind /= -1) then
            next = next + 1

            ! make l the next node of the slice
            nodes(next) = l

            ! and put the corresponding matrix entry into the array `vals`
            vals(next) = val(ind)
        endif
    enddo

end subroutine get_slice_discontiguous




!==========================================================================!
!==== Mutator kernels                                                  ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine graph_leftperm(g,val,p)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val_tmp(:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call g%left_permute(p, edge_p)

    ! If the edges require permutation,
    if (size(edge_p,2) > 0) then
        ! make a temporary array for the matrix values
        allocate( val_tmp(g%capacity) )
        val_tmp = val

        ! Rearrange the values
        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = val_tmp(source:source+num_nodes-1)
        enddo

        ! Deallocate the temporary array
        deallocate(val_tmp)
    endif

    ! Deallocate the array describing the edge permutation
    deallocate(edge_p)

end subroutine graph_leftperm



!--------------------------------------------------------------------------!
subroutine graph_rightperm(g,val,p)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: p(:)
    ! local variables
    integer :: i, source, dest, num_nodes
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val_tmp(:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call g%right_permute(p, edge_p)

    ! If the edges require permutation,
    if (size(edge_p,2) > 0) then
        ! make a temporary array for the matrix values
        allocate( val_tmp(g%capacity) )
        val_tmp = val

        ! Rearrange the values
        do i=1,size(edge_p,2)
            source = edge_p(1,i)
            dest = edge_p(2,i)
            num_nodes = edge_p(3,i)

            val(dest:dest+num_nodes-1) = val_tmp(source:source+num_nodes-1)
        enddo

        ! Deallocate the temporary array
        deallocate(val_tmp)
    endif

    ! Deallocate the array describing the edge permutation
    deallocate(edge_p)

end subroutine graph_rightperm




end module default_sparse_matrix_kernels
