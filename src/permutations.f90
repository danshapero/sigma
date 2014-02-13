!==========================================================================!
!==========================================================================!
module permutations                                                        !
!==========================================================================!
!==========================================================================!
!==== This module contains subroutines for re-ordering graphs to       ====!
!==== accelerate sparse matrix computations, e.g.the Cuthill-McKee     ====!
!==== ordering, multicolor orderings, independent set orderings, etc.  ====!
!==========================================================================!
!==========================================================================!

use graphs
use types, only: linked_list, dynamic_array

implicit none

contains



!--------------------------------------------------------------------------!
subroutine breadth_first_search(g,p)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    integer, intent(out) :: p(:)
    ! local variables
    integer :: i,j,k,num,nbrs(g%max_degree)
    type(linked_list) :: queue

    num = 0

    ! Initialize a queue
    call queue%append(1)

    ! Mark all nodes as unvisited
    p = -1

    ! So long as the queue is not empty,
    do while(queue%length>0)
        num = num+1

        ! Pop the first value i off
        i = queue%get_value(1)
        call queue%delete_entry(1)

        ! Label i according to the order in which it was visited
        p(i) = num

        ! Find all neighbors of node i
        call g%neighbors(i,nbrs)

        ! For each neighbor,
        do k=1,g%max_degree
            j = nbrs(k)
            if (j/=0) then
                ! if that neighbor has not been visited,
                if (p(j)==-1) then
                    ! put it into the queue and label it as such.
                    call queue%append(j)
                    p(j) = 0
                endif
            endif
        enddo
    enddo

end subroutine breadth_first_search



!--------------------------------------------------------------------------!
subroutine greedy_coloring(g,colors)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    integer, intent(out) :: colors(:)
    ! local variables
    integer :: i,j,k,nbrs(g%max_degree),nbrs_colors(g%max_degree+1), &
        & color_totals(g%max_degree+1), used_colors, color, min_occupancy
    type(linked_list) :: queue

    ! Initialize a queue
    call queue%append(1)

    ! Set all colors to a negative value to indicate they have not been
    ! visited or enqueued yet, and set color(1) to 0 to show that it's been
    ! queued up
    colors    = -1
    colors(1) = 0

    ! Set the tally of how many nodes there are of each color
    color_totals = 0
    used_colors  = 0

    ! Keep iterating until the queue is empty
    do while (queue%length>0)
        nbrs_colors = 0

        ! Pop the first value i off the queue
        i = queue%get_value(1)
        call queue%delete_entry(1)

        ! Find all the neighbors of node i
        call g%neighbors(i,nbrs)

        ! Tally up how many neighbors there are of each color
        do k=1,g%max_degree
            j = nbrs(k)
            if (j/=0) then
                color = colors(j)
                if (color>0) then
                    nbrs_colors(color) = nbrs_colors(color)+1
                elseif (color==-1) then
                    !!----------------------------------------------------!!
                    !! When the queue grows large, this operation will be !!
                    !! slow because there are so many pointer redirects   !!
                    !! until we can find the end of the queue on which to !!
                    !! append node j. How can we come up with a better    !!
                    !! queue data structure than a linked list or DA?     !!
                    !! Modify the linked_list type so that we also store  !!
                    !! a tail pointer? Modify the DA type so that we can  !!
                    !! remove from the front as well?                     !!
                    !! Could we prepend it instead and do a depth-first   !!
                    !! rather than breadth-first search?                  !!
                    !!----------------------------------------------------!!
                    call queue%append(j)
                    colors(j) = 0
                endif
            endif
        enddo

        color = 0
        min_occupancy = g%n+1
        ! Find which neighboring color has the fewest nodes assigned to it
        do k=1,used_colors
            if (color_totals(k)>0 .and. color_totals(k)<min_occupancy &
                & .and. nbrs_colors(k)==0) then
                color = k
                min_occupancy = color_totals(k)
            endif
        enddo

        ! If the color chosen is still zero, that means that, of all the
        ! colors we're already using, there is some neighbor that already
        ! has that color. This meas that node i needs to be given a brand
        ! new color we haven't used before.
        if (color==0) then
            color = used_colors+1
            used_colors = used_colors+1
        endif

        colors(i) = color
        color_totals(color) = color_totals(color)+1
    enddo

    call queue%free()

end subroutine greedy_coloring



!--------------------------------------------------------------------------!
subroutine greedy_color_ordering(g,p,ptrs,num_colors)                      !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(graph), intent(in) :: g
    integer, intent(out) :: p(:), ptrs(:), num_colors
    ! local variables
    integer :: i,k,color,added(g%max_degree+1)

    ! Assign colors to all the nodes and put them in the array p
    call greedy_coloring(g,p)

    ! Find the total number of colors used
    num_colors = maxval(p)

    ! Make an array for the starting point in the permutation of each color
    ptrs = 0
    do i=1,g%n
        ptrs(p(i)+1) = ptrs(p(i)+1)+1
    enddo

    ptrs(1) = 1
    do color=1,num_colors
        ptrs(color+1) = ptrs(color+1)+ptrs(color)
    enddo

    ! Build the permutation
    added = 0
    do i=1,g%n
        color = p(i)
        k = added(color)
        p(i) = ptrs(color)+k

        added(color) = added(color)+1
    enddo

end subroutine greedy_color_ordering




end module permutations
