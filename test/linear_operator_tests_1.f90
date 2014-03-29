program linear_operator_tests_1

use sigma

implicit none

    integer :: i,j,k,nn
    class(graph), pointer :: g, h
    real(dp) :: w, z, p
    type(sparse_matrix), pointer :: A, B
    real(dp), allocatable :: x(:), y(:)
    class(linear_operator), pointer :: L
    type(graph_edge_cursor) :: cursor
    integer :: num_blocks, num_returned, edges(2,64)

    call init_seed()

    !---------------------------
    ! Create some random graphs
    !---------------------------
    nn = 256
    p = 4.0/nn
    allocate(ll_graph::g)
    allocate(ll_graph::h)

    call g%init(nn)
    call h%init(nn)

    do i=1,nn
        do j=1,nn
            call random_number(z)
            if (z<p) call g%add_edge(i,j)

            call random_number(z)
            if (z<p) call h%add_edge(i,j)
        enddo
    enddo


    !--------------------------------------
    ! Make random matrices on those graphs
    !--------------------------------------
    allocate(A,B)
    call A%init(nn,nn,'row',g)
    call B%init(nn,nn,'col',h)

    cursor = g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        edges = g%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            call random_number(z)
            call random_number(w)

            call A%set_entry(i,j,dsqrt(-2*log(w))*cos(2*pi*z))
        enddo
    enddo


    cursor = h%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        edges = h%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(2,k)
            j = edges(1,k)

            call random_number(z)
            call random_number(w)

            call B%set_entry(i,j,dsqrt(-2*log(w))*cos(2*pi*z))
        enddo
    enddo


    !--------------------------------------------------------------------
    ! Make a linear operator which is the sum of the two sparse matrices
    !--------------------------------------------------------------------
    L = add_operators(A,B)

    

end program linear_operator_tests_1
