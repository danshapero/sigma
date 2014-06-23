module random_graphs

use types, only: dp
use util, only: init_seed
use graphs
use ellpack_graphs
use cs_graphs

implicit none


contains


!--------------------------------------------------------------------------!
subroutine erdos_renyi(g,nn,p)                                             !
!--------------------------------------------------------------------------!
!     Generate an Erdos-Renyi graph g on n vertices.                       !
!     Any two edges (i,j) in g are connected independently of each other   !
! with probability p.                                                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: nn
    real(dp), intent(in) :: p
    ! local variables
    integer :: i,j
    real(dp) :: z(nn)

    call init_seed()
    call g%init(nn,nn)

    do i=1,nn
        call random_number(z)

        do j=i+1,nn
            if (z(j)<p) then
                call g%add_edge(i,j)
                call g%add_edge(j,i)
            endif
        enddo
    enddo

end subroutine erdos_renyi



!--------------------------------------------------------------------------!
subroutine watts_strogatz(g,nn,k,p)                                        !
!--------------------------------------------------------------------------!
!     Generate a Watts-Strogatz graph g on the n-vertex ring, where each   !
! vertex is connected to the following & preceding k vertices in the ring. !
!     The edges are then independently re-wired with probability p.        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: nn, k
    real(dp), intent(in) :: p
    ! local variables
    integer :: i,j,l,m,n
    real(dp) :: z
    ! variables for graph edge iterator
    type(graph_edge_cursor) :: cursor
    integer :: edges(2,batch_size), num_returned, num_batches
    ! ring graph
    type(ellpack_graph) :: g_ring

    call init_seed()

    !------------------------------
    ! First, generate a ring graph
    call g_ring%init(nn,nn,degree=2*k)

    do i=1,nn
        do m=1,k
            j = mod(i+m-1,nn)+1
            call g_ring%add_edge(i,j)
            call g_ring%add_edge(j,i)
        enddo
    enddo


    !------------------------------------------
    ! Next, rewire the edges of the ring graph
    call g%init(nn)

    cursor = g_ring%make_cursor(0)
    print *, cursor%final,cursor%start
    num_batches = (cursor%final-cursor%start)/batch_size+1
    do n=1,num_batches
        call g_ring%get_edges(edges,cursor,batch_size,num_returned)

        do m=1,num_returned
            ! Get the next edge (i,j) from the ring graph
            i = edges(1,m)
            j = edges(2,m)

            ! Set a vertex l to point to j
            l = j
            if (i/=0 .and. j>i) then
                ! Pick a uniform random number
                call random_number(z)

                ! If it's less than the threshold probability,
                if (z<p) then
                    ! Rewire the edge (i,j); pick a random, new value of
                    ! l so long as it's different from j and we don't
                    ! create a self-loop or multiple edge
                    do while (l==j .or. l==i .or. g%connected(i,l))
                        call random_number(z)
                        l = min(int(z*nn)+1,nn)
                    enddo
                endif
            endif

            ! Whether or not l is still equal to j, add the edge (i,l)
            ! into g
            call g%add_edge(i,l)
            call g%add_edge(l,i)
        enddo
    enddo

end subroutine watts_strogatz



!--------------------------------------------------------------------------!
subroutine barabasi_albert(g,nn,k)                                         !
!--------------------------------------------------------------------------!
!     Generate a scale-free graph g on n vertices using the preferential   !
! attachment algorithm of Barabasi-Albert; P[degree of node = k] = k^(-p). !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(inout) :: g
    integer, intent(in) :: nn, k
    ! local variables
    integer :: i,j,l,d,d_sum
    real(dp) :: z

    call init_seed()
    call g%init(nn)

    ! Generate an initial connected graph on the first k vertices
    do i=1,k-1
        j = i+1

        call g%add_edge(i,j)
        call g%add_edge(j,i)
    enddo

    ! Add new connections to new vertices in succession
    do i=k+1,nn
        ! Add k new connections
        do l=1,k
            d_sum = 2*g%ne

            ! Generate a random number z
            call random_number(z)

            d = 0
            do j=1,i-1
                if (d <= z*d_sum .and. z*d_sum < d+g%degree(j)) then
                    
                    call g%add_edge(i,j)
                    call g%add_edge(j,i)
                endif
            enddo
        enddo
    enddo
    

end subroutine barabasi_albert






end module random_graphs
