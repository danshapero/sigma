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
    slice = 0.0_dp

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

    ! Set the returned nodes to 0
    nodes = 0

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
            slice(next) = val(ind)
        endif
    enddo

end subroutine get_slice_discontiguous




!==========================================================================!
!==== Mutator kernels                                                  ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine set_matrix_value_with_reallocation(g, val, i, j, z)             !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    real(dp), allocatable :: val(:)
    integer, intent(in) :: i, j
    real(dp), intent(in) :: z
    ! local variables
    real(dp), allocatable :: val_temp(:)
    integer :: k, l, m, indx
    integer :: n, num_batches, num_returned, edges(2, batch_size)
    type(graph_edge_cursor) :: cursor
    class(graph), pointer :: h

    allocate(h, mold = g)
    call h%copy(g)

    call g%add_edge(i, j)

    allocate(val_temp(g%capacity))
    indx = g%find_edge(i, j)
    val_temp(indx) = z

    cursor = g%make_cursor()
    num_batches = (cursor%final - cursor%start)/batch_size + 1

    do n = 1, num_batches
        call g%get_edges(edges, cursor, batch_size, num_returned)

        do m = 1, num_returned
            k = edges(1, k)
            l = edges(2, k)
            indx = g%find_edge(k, l)
            val_temp(indx) = val(batch_size * (n - 1) + m)
        enddo
    enddo

    call move_alloc(from = val_temp, to = val)

end subroutine set_matrix_value_with_reallocation



!--------------------------------------------------------------------------!
subroutine graph_leftperm(g,val,p)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: p(:)
    ! local variables
    integer, allocatable :: edge_p(:,:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call g%left_permute(p, edge_p)

    ! Rearrange the array `val` according to the edge permutation
    call rearrange_array_with_compressed_permutation(val, edge_p)

    ! Deallocate the array describing the edge permutation
    deallocate(edge_p)

end subroutine graph_leftperm



!--------------------------------------------------------------------------!
subroutine graph_rightperm(g, val, p)                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: p(:)
    ! local variables
    integer, allocatable :: edge_p(:,:)

    ! Permute the graph and get an array edge_p describing the permutation
    ! of the edges
    call g%right_permute(p, edge_p)

    ! Rearrange the array `val` according to the edge permutation
    call rearrange_array_with_compressed_permutation(val, edge_p)

    ! Deallocate the permutation array
    deallocate(edge_p)

end subroutine graph_rightperm



!--------------------------------------------------------------------------!
subroutine assemble_matrix(g, val)                                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph) :: g
    real(dp), allocatable :: val(:)
    ! local variables
    integer :: k, source, dest, num
    integer, allocatable :: edge_p(:,:)
    real(dp), allocatable :: val_temp(:)

    call g%compress(edge_p)

    if (size(edge_p, 2) /= 0) then
        ! Make a temporary array for the matrix entries
        allocate( val_temp(g%capacity) )

        do k = 1, size(edge_p, 2)
            source = edge_p(1, k)
            dest   = edge_p(2, k)
            num    = edge_p(3, k)

            val_temp(dest : dest + num - 1) = val(source : source + num - 1)
        enddo

        call move_alloc(from = val_temp, to = val)
    endif

end subroutine assemble_matrix




!==========================================================================!
!==== Auxiliary routines                                               ====!
!==========================================================================!

!--------------------------------------------------------------------------!
subroutine rearrange_array_with_compressed_permutation(val, p)             !
!--------------------------------------------------------------------------!
!     This is a commonly used routine. A compressed permutation `p` is an  !
! array of dimension `p(3, m)`; it specifies how an array `val` of length  !
! `n` is to be permuted according to the following rule:                   !
!     p(1, k) = source index in `val` of block `k`                         !
!     p(2, k) = destination index in the permuted array of block `k`       !
!     p(3, k) = size of block `k`.                                         !
! In this fashion, a permutation can be represented with much less space   !
! when we are moving around large contiguous chunks of the array. For      !
! example, the compressed permutation                                      !
!     [ 1, n/2 + 1, n/2;                                                   !
!       n/2 + 1, 1, n/2 ]                                                  !
! represents swapping the first and second halves of an array with only 6  !
! integers, no matter what the array size `n` is.                          !
!     In general, so long as the average block size to move is greater     !
! 3, this routine is faster than the usual approach.                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    real(dp), intent(inout) :: val(:)
    integer, intent(in) :: p(:,:)
    ! local variables
    integer :: k, source, dest, num
    real(dp) :: val_temp(size(val))

    if (size(p, 2) /= 0) then
        val_temp = val

        do k = 1, size(p, 2)
            source = p(1, k)
            dest   = p(2, k)
            num    = p(3, k)

            val(dest : dest + num - 1) = val_temp(source : source + num - 1)
        enddo
    endif

end subroutine rearrange_array_with_compressed_permutation



end module default_sparse_matrix_kernels
