program linear_operator_tests_1

use sigma

implicit none

    integer :: i,j,k,n,nn
    class(graph), pointer :: g, h
    real(dp) :: p, q, r
    type(sparse_matrix), pointer :: A, B
    real(dp), allocatable :: w(:), x(:), y(:), z(:)
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
            call random_number(q)
            if (q<p) call g%add_edge(i,j)

            call random_number(q)
            if (q<p) call h%add_edge(i,j)
        enddo
    enddo


    !--------------------------------------
    ! Make random matrices on those graphs
    !--------------------------------------
    allocate(A,B)
    call A%init(nn,nn,'row',g)

    cursor = g%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        edges = g%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            call random_number(q)
            call random_number(r)

            call A%set_value(i,j,dsqrt(-2*log(r))*cos(2*pi*q))
        enddo
    enddo

    call B%init(nn,nn,'col',h)
    cursor = h%make_cursor(0)
    num_blocks = (cursor%final-cursor%start)/64+1

    do n=1,num_blocks
        edges = h%get_edges(cursor,64,num_returned)

        do k=1,num_returned
            i = edges(2,k)
            j = edges(1,k)

            call random_number(q)
            call random_number(r)

            call B%set_value(i,j,dsqrt(-2*log(r))*cos(2*pi*q))
        enddo
    enddo


    !--------------------------------------------------------------------
    ! Make a linear operator which is the sum of the two sparse matrices
    !--------------------------------------------------------------------
    L => add_operators(A,B)



    !--------------------------
    ! Multiply L by some stuff
    !--------------------------
    allocate(w(nn),x(nn),y(nn),z(nn))
    x = 1.0_dp

    z = 0.0_dp
    call A%matvec_add(x,z)
    call B%matvec_add(x,z)
    if (maxval(dabs(z))==0.0_dp) then
        print *, 'Something went way wrong, made random matrices A, B'
        print *, 'and (A+B)*[1,...,1] = 0. Terminating.'
        call exit(1)
    endif

    y = 0.0_dp
    call L%matvec_add(x,y)
    if (maxval(dabs(y))==0.0_dp) then
        print *, 'Setting L = A+B and multiplying y = L*[1,..,1] failed,'
        print *, 'got y = 0. Terminating.'
        call exit(1)
    endif

    r = maxval(dabs(y-z))
    if (r>1.0e-14) then
        print *, 'Setting L = A+B and multiplying y = L*[1,..,1] failed;'
        print *, 'computed z = (A*x)+(B*x) but ||y-z|| = ',r
        print *, 'Terminating.'
        call exit(1)
    endif

    nullify(L)


    !--------------------------------------------------------------------
    ! Make a linear operator which is the product of two sparse matrices
    !--------------------------------------------------------------------
    L => multiply_operators(A,B)

    x = 1.0_dp
    call A%matvec(x,w)
    call B%matvec(w,z)
    if (maxval(dabs(z))==0.0_dp) then
        print *, 'Something went way wrong, made random matrices A, B'
        print *, 'and A*B*[1,..,1] = 0. Terminating.'
        call exit(1)
    endif

    y = 0.0_dp
    call L%matvec(x,y)
    if (maxval(dabs(y))==0.0_dp) then
        print *, 'Setting L = A*B and multiplying y = L*[1,..,1] failed,'
        print *, 'got y = 0. Terminating.'
        call exit(1)
    endif

    r = maxval(dabs(y-z))
    if (r>1.0e-14) then
        print *, 'Setting L = A*B and multiplying y = L*[1,..,1] failed;'
        print *, 'computed z = A*(B*x) but ||y-z|| = ',r
        print *, 'Terminating.'
        call exit(1)
    endif

end program linear_operator_tests_1
