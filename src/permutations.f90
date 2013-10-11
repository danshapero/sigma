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




end module permutations
