module permutation_mod

    use sparse_matrix_mod

    implicit none

contains



!--------------------------------------------------------------------------!
subroutine bfs(A,p)                                                        !
!--------------------------------------------------------------------------!
! Compute the breadth-first search ordering for the matrix A. This order   !
! is useful for problems where one intends to use an incomplete Cholesky   !
! preconditioner, as the fill-in tends to be dramatically reduced.         !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: p(A%nrow)
    ! local variables
    integer :: i,j,col,next,q(A%nrow),nbrs(A%max_degree)

    ! q is an auxiliary array, which is the inverse of the permutation by
    ! the time we're done computing it. Using this array makes for a faster
    ! lookup of whether or not a node has been added to the permutation
    ! than searching through it every time.
    q = 0
    p = 0
    q(1) = 1
    p(1) = 1
    next = 1

    do i=1,A%nrow
        nbrs = A%get_neighbors( q(i) )
        do j=1,A%max_degree
            if (nbrs(j)/=0) then
                if ( p(nbrs(j))==0 .and. nbrs(j)/=i ) then
                    next = next+1
                    q(next) = nbrs(j)
                    p( nbrs(j) ) = next
                endif
            endif
        enddo

    enddo

    ! Reverse the order of the permutation
    do i=1,A%nrow
        p(i) = A%nrow-p(i)+1
    enddo

end subroutine bfs



!--------------------------------------------------------------------------!
subroutine greedy_multicolor(A,p,maxcolor)                                 !
!--------------------------------------------------------------------------!
! Compute a greedy multi-color ordering for the matrix A.                  !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: maxcolor,p(A%nrow)
    ! local variables
    integer :: i,j,clr,colors(A%nrow),nbrs(A%max_degree), &
        & num_colors(A%max_degree+1)

    nbrs = 0
    colors = 0
    maxcolor = 1
    p = 0
    p(1) = 1

    do i=1,A%nrow
        num_colors = 0
        nbrs = A%get_neighbors(i)
        do j=1,A%max_degree
            if (nbrs(j)/=0) then
                clr = colors(nbrs(j))
                if (clr/=0) then
                    num_colors(clr) = num_colors(clr)+1
                endif
            endif
        enddo
        do clr=A%max_degree+1,1,-1
            if (num_colors(clr)==0) colors(i) = clr
        enddo
    enddo
    maxcolor = maxval(colors)

    num_colors = 0
    do i=1,A%nrow
        clr = colors(i)
        num_colors(clr) = num_colors(clr)+1
    enddo

    do i=maxcolor,2,-1
        num_colors(i) = sum(num_colors(1:i-1))+1
    enddo
    num_colors(1) = 1

    do i=1,A%nrow
        clr = colors(i)
        p(i) = num_colors(clr)
        num_colors(clr) = num_colors(clr)+1
    enddo


end subroutine greedy_multicolor



!--------------------------------------------------------------------------!
subroutine block_multicolor(A,p,blocks)                                    !
!--------------------------------------------------------------------------!
    implicit none
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    integer, intent(out) :: p(:),blocks(:)
    ! local variables
    integer :: i,j,k,l,m,n,colors(A%nrow),nbrs(A%max_degree), &
        & num_colors(A%max_degree+1),next(A%nrow),num_nodes_colored, &
        & head,tail,maxcolor

    nbrs = 0
    colors = 0
    maxcolor = 1
    p = 0
    next = 0
    p(1) = 1
    head = 1
    tail = 1
    num_nodes_colored = 0

    do while(num_nodes_colored < A%nrow)
        nbrs = A%get_neighbors(head)
        do i=1,A%max_degree
            n = nbrs(i)
            if ( n/=0 .and. n/=head ) then
                if ( next(n)/=0 ) then
                    next(tail) = n
                    tail = n
                endif
            endif
        enddo
        head = next(head)
        num_nodes_colored = num_nodes_colored+1
    enddo


end subroutine block_multicolor




end module permutation_mod
