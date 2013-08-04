program matrix_tests

use csr_matrices
use sparse_matrices
use graphs

implicit none

    class(graph), allocatable :: g
    class(sparse_matrix), allocatable :: A
    integer :: i,j
    integer, allocatable :: edges(:,:)

    allocate(edges(2,24))

    edges(:,1) = [1, 2]
    edges(:,2) = [1, 3]
    edges(:,3) = [1, 4]
    edges(:,4) = [1, 5]
    edges(:,5) = [1, 6]
    edges(:,6) = [1, 7]
    edges(:,7) = [2, 3]
    edges(:,8) = [3, 4]
    edges(:,9) = [4, 5]
    edges(:,10) = [5, 6]
    edges(:,11) = [6, 7]
    edges(:,12) = [7, 2]

    do i=1,12
        edges(1,i+12) = edges(2,i)
        edges(2,i+12) = edges(1,i)
    enddo

    allocate(cs_graph::g)
    call g%init(7,7,edges)

    allocate(csr_matrix::A)
    call A%assemble(g)

end program matrix_tests
