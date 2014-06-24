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
    integer :: num_batches, num_returned, edges(2,batch_size)

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

    cursor = g%make_cursor()
    num_batches = (cursor%final-cursor%start)/batch_size+1

    do n=1,num_batches
        call g%get_edges(edges,cursor,batch_size,num_returned)

        do k=1,num_returned
            i = edges(1,k)
            j = edges(2,k)

            call random_number(q)
            call random_number(r)

            call A%set_value(i,j,dsqrt(-2*log(r))*cos(2*pi*q))
        enddo
    enddo

    call B%init(nn,nn,'col',h)
    cursor = h%make_cursor()
    num_batches = (cursor%final-cursor%start)/batch_size+1

    do n=1,num_batches
        call h%get_edges(edges,cursor,batch_size,num_returned)

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
    L = A+B

    ! Check that getting the entries of an operator sum works properly
    do k=1,nn
        call random_number(q)
        call random_number(r)

        i = min(int(q*nn)+1,nn)
        j = min(int(r*nn)+1,nn)

        q = dabs(L%get_value(i,j) - (A%get_value(i,j)+B%get_value(i,j)))

        if (q>1.0e-12) then
            print *, 'Getting entry',i,j,'from operator sum failed;'
            print *, 'A+B = ',A%get_value(i,j)+B%get_value(i,j)
            print *, 'L = ',L%get_value(i,j)
            call exit(1)
        endif
    enddo

    ! Check that operator-vector multiplication works properly if
    ! an operator is defined as the sum of two operators
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
    L = A*B

    x = 1.0_dp
    call B%matvec(x,w)
    call A%matvec(w,z)
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

    nullify(L)



    !----------------------------------------------------------------------!
    ! Make a linear operator which is the adjoint of another operator      !
    !----------------------------------------------------------------------!
    L = adjoint(A)

    x = 1.0_dp
    call A%matvec(x,y,.true.)
    call L%matvec(x,z)

    r = maxval(dabs(y-z))
    if (r>1.0e-12) then
        print *, 'Constructing adjoint L = A* of operator A failed;'
        print *, 'For x = [1,..,1], || Lx - A*x || =',r
        print *, 'Terminating.'
        call exit(1)
    endif

    nullify(L)



    !----------------------------------------------------------------------!
    ! Form the product A*A of an operator                                  !
    !----------------------------------------------------------------------!
    L = adjoint(A)*A

    x = 1.0_dp
    call L%matvec(x,y)
    call A%matvec(x,w)
    call A%matvec(w,z,.true.)

    r = maxval(dabs(y-z))
    if (r>1.0e-12) then
        print *, 'Constructor operator L = A*A failed;'
        print *, 'For x = [1,..,1], || Lx - A*Ax || =',r
        print *, 'Terminating.'
        call exit(1)
    endif

end program linear_operator_tests_1
