program matrix_tests_3

use sigma

implicit none

    type(sparse_matrix) :: A
    class(graph), pointer :: g

    integer :: i,j,nn


    nn = 64

    !----------------------------------------------------------------------!
    ! Construct graph for connectivity structure of sparse matrix          !
    !----------------------------------------------------------------------!
    allocate(cs_graph::g)
    call g%init(nn,degree=3)

    do i=1,nn
        call g%add_edge(i,i)

        j = mod(i,nn)+1
        call g%add_edge(i,j)

        j = mod(i+nn-2,nn)+1
        call g%add_edge(i,j)
    enddo


    !----------------------------------------------------------------------!
    ! Make a matrix with g as its connectivity structure                   !
    !----------------------------------------------------------------------!
    call A%init(nn,nn,'row',g)
    call A%zero()

    do i=1,nn
        call A%set_value(i,i,5.0_dp)

        j = mod(i,nn)+1
        call A%set_value(i,j,-1.0_dp)
        call A%set_value(j,i,-1.0_dp)
    enddo

    if (A%g%max_degree/=3) then
        print *, 'Max degree of A%g should be 3 after initializing and'
        print *, 'filling circulant matrix.'
        print *, 'Value found:',A%g%max_degree
        print *, 'Terminating'
        call exit(1)
    endif

    if (A%g%connected(1,5) .or. A%g%connected(1,7)) then
        print *, 'A%g should not have 1,5 or 1,7 connected yet.'
        print *, 'Terminating.'
        call exit(1)
    endif



    !----------------------------------------------------------------------!
    ! Add some entries to A that weren't present in g                      !
    !----------------------------------------------------------------------!
    do i=1,nn
        j = mod(i+3,nn)+1

        call A%set_value(i,j,-1.0_dp)
        call A%set_value(j,i,-1.0_dp)
    enddo

    if (A%g%max_degree/=5) then
        print *, 'Max degree of A%g should be 5 after adding 2 matrix '
        print *, 'entries to every row where they were not already pre-'
        print *, 'allocated. Value found:',A%g%max_degree
        print *, 'Terminating.'
        call exit(1)
    endif

    if (A%nnz/=5*nn) then
        print *, 'Number of non-zero entries of A after adding 2 matrix'
        print *, 'entries to every row where they were not already pre-'
        print *, 'allocated should be',5*nn
        print *, 'Value found:',A%nnz
        print *, 'Terminating.'
        call exit(1)
    endif


    do i=1,nn
        j = mod(i+5,nn)+1

        call A%set_value(i,j,-1.0_dp)
        call A%set_Value(j,i,-1.0_dp)
    enddo


end program matrix_tests_3
