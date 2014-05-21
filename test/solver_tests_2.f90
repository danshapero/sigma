program solver_tests_2

use sigma

implicit none

    ! graph and matrix
    class(graph), pointer :: g
    type(sparse_matrix) :: A
    ! vectors
    real(dp), allocatable :: x(:), xc(:), b(:), r(:), w(:), q(:)
    ! LDU solver
    class(sparse_ldu_solver), pointer :: ildu
    ! indices and miscellaneous
    integer :: i,j,k,d,nn,iter
    integer, allocatable :: neighbors(:)
    real(dp) :: p, z, misfit
    ! command-line arguments
    character(len=16) :: arg
    logical :: verbose


    nn = 128
    p = log(1.0_dp*nn)/nn
    call init_seed()

    !----------------------------------------------------------------------!
    ! Get command line arguments to see if we're running in verbose mode   !
    !----------------------------------------------------------------------!
    verbose = .false.
    call getarg(1,arg)
    if (trim(arg)=="-v" .or. trim(arg)=="--verbose") then
        verbose = .true.
    endif


    !----------------------------------------------------------------------!
    ! Create a random graph                                                !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)
    call g%init(nn)
    do i=1,nn
        call g%add_edge(i,i)

        do j=i+1,nn
            call random_number(z)

            if (z<p) then
                call g%add_edge(i,j)
                call g%add_edge(j,i)
            endif
        enddo
    enddo

    call g%compress()

    if (verbose) then
        print *, 'o Done creating random connectivity graph G.'
        print *, '    Nodes:',nn
        print *, '    Edges:',g%ne/2
        print *, '    Maximum degree:',g%max_degree
    endif


    !----------------------------------------------------------------------!
    ! Make the Laplacian of the graph                                      !
    !----------------------------------------------------------------------!
    call A%init(nn,nn,'row',g)

    allocate(neighbors(g%max_degree))

    do i=1,nn
        call A%set_value(i,i,1.0_dp)

        call g%get_neighbors(neighbors,i)

        d = g%degree(i)
        do k=1,d
            j = neighbors(k)

            if (j/=i) then
                call A%set_value(i,j,-1.0_dp)
                call A%add_value(i,i,+1.0_dp)
            endif
        enddo
    enddo

    if (verbose) print *, 'o Done creating A = I + L(G).'


    !----------------------------------------------------------------------!
    ! Build a direct solver for the system                                 !
    !----------------------------------------------------------------------!
    allocate(ildu)
    call ildu%setup(A)

    if (verbose) then
        print *, 'o Done constructing incomplete LDU'
        print *, '    factorization of A.'
    endif


    !----------------------------------------------------------------------!
    ! Solve a linear system with the direct solver                         !
    !----------------------------------------------------------------------!
    allocate(x(nn), xc(nn), b(nn), r(nn), q(nn), w(nn))
    x = 0.0_dp
    call random_number(xc)
    call A%matvec(xc,b)

    if (verbose) then
        print *, 'o Creating random right-hand side b with'
        print *, '    which to solve system A*x = b.'
    endif

    r = b
    do iter=1,nn
        call ildu%solve(A,q,r)
        x = x+q
        call A%matvec(x,w)
        r = b-w
    enddo

    misfit = maxval(dabs(x-xc))

    if (misfit > 1.0e-15) then
        print *, 'Stationary iterative method with incomplete'
        print *, 'LDU decomposition did not converge.'
        print *, 'Terminating.'
        call exit(1)
    endif

    if (verbose) then
        print *, 'o Done solving linear system.'
        print *, '    Error:',misfit
    endif
    

end program solver_tests_2
