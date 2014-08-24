!--------------------------------------------------------------------------!
program linear_operator_test_algebra                                       !
!--------------------------------------------------------------------------!
! This program tests algebraic operations on linear operators:             !
!     o sums                                                               !
!     o products                                                           !
!     o adjoints                                                           !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs
    class(graph_interface), pointer :: g, h

    ! sparse matrices
    type(csr_matrix), target :: A
    type(csc_matrix), target :: B

    ! linear operators
    class(linear_operator), pointer :: L

    ! vectors
    real(dp), dimension(:), allocatable :: w, x, y, z

    ! integer indices
    integer :: i, j, k, n, nn

    ! random variables
    real(dp) :: p, q, r

    ! variables for iterating over matrix entries
    type(graph_edge_cursor) :: cursor
    integer :: num_batches, num_returned, edges(2, batch_size)

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose



    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!
    verbose = .false.
    call getarg(1,arg)
    select case(trim(arg))
        case("-v")
            verbose = .true.
        case("-V")
            verbose = .true.
        case("--verbose")
            verbose = .true.
    end select



    !----------------------------------------------------------------------!
    ! Set the matrix size and initialize a random seed                     !
    !----------------------------------------------------------------------!

    nn = 64
    p = log(1.0_dp * nn) / log(2.0_dp) / nn

    call init_seed()



    !----------------------------------------------------------------------!
    ! Make random graphs                                                   !
    !----------------------------------------------------------------------!

    allocate(ll_graph :: g)
    allocate(ll_graph :: h)

    call g%init(nn)
    call h%init(nn)

    do i = 1, nn
        do j = 1, nn
            call random_number(q)
            if (q < p) call g%add_edge(i, j)

            call random_number(q)
            if (q < p) call h%add_edge(j, i)
        enddo
    enddo

    if (verbose) then
        print *, 'Done building random graphs.'
        print *, '    Number of vertices:', nn
        print *, '    Number of edges:   ', g%ne, h%ne
        print *, '    Max vertex degree: ', g%max_degree(), h%max_degree()
    endif



    !----------------------------------------------------------------------!
    ! Make random matrices on those graphs                                 !
    !----------------------------------------------------------------------!

    call A%init(nn, nn, g)
    call B%init(nn, nn, h)


    cursor = g%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call g%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            call random_number(q)
            call A%set_value(i, j, 2 * q - 1)
        enddo
    enddo


    cursor = h%make_cursor()
    num_batches = (cursor%final - cursor%start) / batch_size + 1

    do n = 1, num_batches
        call h%get_edges(edges, cursor, batch_size, num_returned)

        do k = 1, num_returned
            i = edges(1, k)
            j = edges(2, k)

            call random_number(q)

            call B%set_value(i, j, 2 * q - 1)
        enddo
    enddo

    if (verbose) print *, 'Done building random matrices'



    !----------------------------------------------------------------------!
    ! Make a linear operator as the sum of `A` and `B`                     !
    !----------------------------------------------------------------------!

    L = A + B

    if (verbose) print *, 'Testing operator sum'


    ! Check that getting the entries of an operator sum works
    do i = 1, nn
        do j = 1, nn
            q = dabs( L%get_value(i, j) &
                    & - A%get_value(i, j) - B%get_value(i, j) )
            if (q > 1.0e-14) then
                print *, 'Getting entry', i, j, 'from operator sum failed.'
            call exit(1)
            endif
        enddo
    enddo


    ! Check that multiplying an operator sum by a vector works
    allocate( w(nn), x(nn), y(nn), z(nn) )
    x = 1.0_dp

    z = 0.0_dp
    call A%matvec_add(x, z)
    call B%matvec_add(x, z)
    if (maxval(dabs(z)) == 0.0_dp) then
        print *, 'Something went way wrong, made random matrices A, B'
        print *, 'and (A+B)*[1,...,1] = 0. Terminating.'
        call exit(1)
    endif

    ! First, check that the result is non-zero
    y = 0.0_dp
    call L%matvec_add(x, y)
    if (maxval(dabs(y)) == 0.0_dp) then
        print *, 'Setting L = A+B and multiplying y = L*[1,..,1] failed,'
        print *, 'got y = 0.'
        call exit(1)
    endif

    ! Then check that it's actually right
    r = maxval(dabs(y - z))
    if (r > 1.0e-14) then
        print *, 'Setting L = A+B and multiplying y = L*[1,..,1] failed;'
        print *, 'computed z = (A*x)+(B*x) but ||y-z|| = ', r
        call exit(1)
    endif

    nullify(L)



    !----------------------------------------------------------------------!
    ! Make a linear operator as the product of `A` and `B`                 !
    !----------------------------------------------------------------------!

    L = A * B

    if (verbose) print *, 'Testing operator product'

    x = 1.0_dp
    call B%matvec(x, w)
    call A%matvec(w, z)

    if (maxval(dabs(z)) == 0.0_dp) then
        print *, 'Something went way wrong, made random matrices A, B'
        print *, 'and A*B*[1,..,1] = 0. Terminating.'
        call exit(1)
    endif


    ! First, check that the result is non-zero
    y = 0.0_dp
    call L%matvec(x, y)

    if (maxval(dabs(y)) == 0.0_dp) then
        print *, 'Setting L = A*B and multiplying y = L*[1,..,1] failed,'
        print *, 'got y = 0. Terminating.'
        call exit(1)
    endif

    ! Then check that it's right
    r = maxval(dabs(y - z))
    if (r > 1.0e-14) then
        print *, 'Setting L = A*B and multiplying y = L*[1,..,1] failed;'
        print *, 'computed z = A*(B*x) but ||y-z|| = ',r
        print *, 'Terminating.'
        call exit(1)
    endif

    nullify(L)



    !----------------------------------------------------------------------!
    ! Make a linear operator as the transpose of `A`                       !
    !----------------------------------------------------------------------!

    L = adjoint(A)

    if (verbose) print *, 'Testing operator adjoint'

    x = 1.0_dp
    call L%matvec(x, y)
    call A%matvec_t(x, z)

    r = maxval(dabs(y - z))
    if (r > 1.0e-12) then
        print *, 'Constructing adjoint L = A* of operator A failed;'
        print *, 'For x = [1,..,1], || Lx - A*x || =', r
        print *, 'Terminating.'
        call exit(1)
    endif

    nullify(L)



    !----------------------------------------------------------------------!
    ! Form the product `adjoint(A) * A` of the matrix `A`                  !
    !----------------------------------------------------------------------!

    L = adjoint(A) * A

    if (verbose) print *, 'Testing adjoint and product'

    x = 1.0_dp
    call L%matvec(x, y)
    call A%matvec(x, w)
    call A%matvec_t(w, z)

    r = maxval(dabs(y - z))
    if (r > 1.0e-12) then
        print *, 'Constructor operator L = A*A failed;'
        print *, 'For x = [1,..,1], || Lx - A*Ax || =',r
        print *, 'Terminating.'
        call exit(1)
    endif


    call A%destroy()
    call B%destroy()
    deallocate(L)
    

end program linear_operator_test_algebra
