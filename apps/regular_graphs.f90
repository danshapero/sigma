module regular_graphs

use graphs

implicit none


contains


!--------------------------------------------------------------------------!
subroutine torus(g,nx,ny)                                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: nx, ny
    ! local variables
    integer :: i,j,x,y

    call g%init(nx*ny,degree=4)

    do x=1,nx
        do y=1,ny
            i = ny*(x-1)+y

            j = ny*(x-1)+mod(y,ny)+1
            call g%add_edge(i,j)
            call g%add_edge(j,i)

            j = ny*mod(x,nx)+y
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo

end subroutine torus



!--------------------------------------------------------------------------!
subroutine petersen(g,n,k)                                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: n, k
    ! local variables
    integer :: i,j,m

    call g%init(2*n,degree=3)

    do i=1,n
        j = mod(i-1,n)+1
        call g%add_edge(i,j)
        call g%add_edge(j,i)

        j = i+n
        call g%add_edge(i,j)
        call g%add_edge(j,i)

        j = mod(i+k-1,n)+1
        call g%add_edge(i+n,j+n)
        call g%add_edge(j+n,i+n)
    enddo

end subroutine petersen



!--------------------------------------------------------------------------!
subroutine flower_snark(g,n)                                               !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: n
    ! local variables
    integer :: i,j,k,l

    call g%init(4*n,degree=3)

    ! Build n copies of the star graph on 4 vertices; for each k<=n,
    ! call the vertices A_i, B_i, C_i and D_i, with A_i the center of
    ! the star
    do k=1,n
        i = 4*(k-1)+1

        do l=1,3
            j = 4*(k-1)+1+l
            call g%add_edge(i,j)
            call g%add_edge(j,i)
        enddo
    enddo

    ! Construct the cycle (B_1,...,B_n)
    do k=1,n
        i = 4*(k-1)+2
        j = mod(4*k+2-1,4*n)+1

        call g%add_edge(i,j)
        call g%add_edge(j,i)
    enddo

    ! Construct the path (C_1,...,C_n)
    do k=1,n-1
        i = 4*(k-1)+3
        j = 4*k+3

        call g%add_edge(i,j)
        call g%add_edge(j,i)
    enddo

    ! Add the edge (C_n,D_1)
    call g%add_edge(4*n-1,4)
    call g%add_edge(4,4*n-1)

    ! Construct the path (D_1,...,D_n)
    do k=1,n-1
        i = 4*k
        j = 4*(k+1)

        call g%add_edge(i,j)
        call g%add_edge(j,i)
    enddo

    ! Add the edge (D_n,C_1)
    call g%add_edge(4*n,3)
    call g%add_edge(3,4*n)

end subroutine flower_snark



!--------------------------------------------------------------------------!
subroutine hypercube(g,n)                                                  !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: n
    ! local variables
    integer :: i,j,k,mask,nn

    nn = 2**n
    call g%init(nn, degree=n)

    do i=1,nn
        do k=1,n
            ! First, make a `mask`, a number 
            !     0...010...0
            ! where the 1 is in the k-th bit.
            mask = 2**(k-1)

            ! If we take the `exclusive or` operation with i-1 and the mask,
            ! we get a new number which is the same as i-1 only with the
            ! k-th bit flipped. Add 1 to this and we get the vertex
            j = ieor(i-1,mask)+1

            call g%add_edge(i,j)
        enddo
    enddo

end subroutine hypercube



end module regular_graphs
