!--------------------------------------------------------------------------!
module permutation_mod                                                     !
!--------------------------------------------------------------------------!
! This module contains subroutines for permuting the entries of sparse     !
! matrices and various sundry matrix reordering procedures, such as the    !
! reverse Cuthill-McKee, multicolor and independent set orderings.         !
!--------------------------------------------------------------------------!
! >> Procedure index                                                       !
! permute_matrix : permute the entries of a sparse matrix                  !
! bfs : breadth-first search through a sparse matrix; the BFS ordering is  !
!       the order in which the nodes were visited                          !
!--------------------------------------------------------------------------!

    use matrix_mod
    use mesh_mod

    implicit none

contains



!--------------------------------------------------------------------------!
subroutine permute_matrix(A,B,permutation)                                 !
!--------------------------------------------------------------------------!
! Perform a symmetric permutation of the entries of a sparse_matrix        !
! Input: a sparse_matrix A, to be permuted                                 !
!        an integer array permutation of length equal to the dimension of  !
!        A; permutation(i) is the destination of node i in the new order   !
! Output: a sparse_matrix B, containing the permuted A                     !
!--------------------------------------------------------------------------!
    implicit none
    ! local variables
    type (sparse_matrix), intent(in) :: A
    type (sparse_matrix), intent(out) :: B
    integer, dimension( A%nn ), intent(in) :: permutation
    ! input/output variables
    integer :: i,j,astart,bstart
    integer, dimension( A%nn ) :: q

    do i=1,A%nn
        q( permutation(i) ) = i
    enddo

    if (.not.allocated(B%ia)) B = init_matrix( A%nn,A%nnz )
    B%ia(1) = 1
    do i=1,A%nn
        bstart = B%ia(i)
        astart = A%ia(q(i))
        B%ia(i+1) = bstart+A%ia(q(i)+1)-A%ia(q(i))
        do j=0,B%ia(i+1)-B%ia(i)-1
            B%ja(bstart+j) = permutation( A%ja(astart+j) )
            B%val(bstart+j) = A%val(astart+j)
        enddo
    enddo

    call sort_ja(B)

end subroutine permute_matrix



!--------------------------------------------------------------------------!
subroutine bfs(A,permutation)                                              !
!--------------------------------------------------------------------------!
! Compute the breadth-first search ordering for the matrix A. This order   !
! is useful for problems where one intends to use an incomplete Cholesky   !
! preconditioner, as the fill-in tends to be dramatically reduced.         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (sparse_matrix), intent(in) :: A
    integer, dimension( A%nn ), intent(out) :: permutation
    ! local variables
    integer :: i,j,col,next
    integer, dimension( A%nn ) :: q

    ! q is an auxiliary array, which is the inverse of the permutation by
    ! the time we're done computing it. Using this array makes for a faster
    ! lookup of whether or not a node has been added to the permutation
    ! than searching through it every time.
    q = 0
    permutation = 0
    q(1) = 1
    permutation(1) = 1
    next = 1

    do i=1,A%nn
        do j=A%ia( q(i) ),A%ia( q(i)+1 )-1
            col = A%ja(j)
            if (permutation(col) == 0) then
                next = next+1
                q(next) = col
                permutation(col) = next
            endif
        enddo
    enddo

    ! Reverse the order of the permutation. Can't remember why I do this
    do i=1,A%nn
        permutation(i) = A%nn-permutation(i)+1
    enddo

end subroutine bfs



!--------------------------------------------------------------------------!
subroutine greedy_multicolor(A,permutation,maxcolor)                       !
!--------------------------------------------------------------------------!
! Compute a greedy multi-color ordering for the matrix A.                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    type (sparse_matrix), intent(in) :: A
    integer, dimension( A%nn ), intent(out) :: permutation
    integer, intent(out) :: maxcolor
    ! local variables
    integer :: i,j
    integer, dimension(A%nn) :: color


end subroutine greedy_multicolor




end module permutation_mod
