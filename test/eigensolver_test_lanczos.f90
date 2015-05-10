!--------------------------------------------------------------------------!
program eigensolver_test_lanczos                                           !
!--------------------------------------------------------------------------!
! This program tests the correctness of the implementation of the Lanczos  !
! algorithm and approximate eigensolvers.                                  !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    type(ll_graph) :: g

    ! sparse matrices
    type(sparse_matrix) :: A

    ! vectors
    real(dp), allocatable :: x(:), y(:)

    ! approximation error
    real(dp) :: err

    ! variables for getting graph neighbors
    integer :: d
    integer, allocatable :: nodes(:)

    ! Lanczos vectors and tridiagonal matrix
    real(dp), allocatable :: V(:,:), alpha(:), beta(:), Q(:,:)

    ! random numbers
    real(dp) :: w, z, p

    ! ints
    integer :: i, j, k, nn, nq

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
    nn = 128
    p = log(1.0_dp * nn) / log(2.0_dp) / nn
    call init_seed()


    !----------------------------------------------------------------------!
    ! Make a random Erdos-Renyi graph                                      !
    !----------------------------------------------------------------------!
    call g%init(nn)

    do i = 1, nn
        call g%add_edge(i, i)

        do j = i + 1, nn
            call random_number(z)

            if (z < p) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo
    enddo

    d = g%get_max_degree()
    allocate(nodes(d))

    if (verbose) then
        print *, "o Done generating Erdos-Renyi graph."
        print *, "    Number of vertices: ", nn
        print *, "    Number of edges:    ", g%get_num_edges()
        print *, "    Maximum degree:     ", d
    endif


    !----------------------------------------------------------------------!
    ! Create a CSR matrix representing the Laplacian of `g`                !
    !----------------------------------------------------------------------!
    call A%set_matrix_type("csr")
    call A%set_dimensions(nn, nn)
    call A%copy_graph(g)
    call A%zero()

    do i = 1, nn
        d = g%get_degree(i)
        call g%get_neighbors(nodes, i)

        do k = 1, d
            j = nodes(k)
            call A%add_value(i, j, -1.0_dp)
            call A%add_value(i, i, +1.0_dp)
        enddo
    enddo

    if (verbose) then
        print *, "o Done creating graph Laplacian."
    endif


    !----------------------------------------------------------------------!
    ! Perform a few steps of the Lanczos process                           !
    !----------------------------------------------------------------------!
    nq = int(dsqrt(1.0_dp * nn))

    allocate(alpha(nq), beta(nq), V(nn, nq), Q(nq, nq), x(nn), y(nn))
    call lanczos(A, alpha, beta, V)

    if (verbose) then
        print *, "o Done executing Lanczos algorithm."
    endif

    do i = 2, nq - 1
        call A%matvec(V(:, i), x)
        y = alpha(i) * V(:, i) + beta(i) * V(:, i-1) + beta(i+1) * V(:, i+1)

        err = dsqrt(sum((y - x) * (y - x)) / sum(x * x))

        if (err > 1.0e-14) then
            print *, "Computing Lanczos vector failed!"
            print *, "Relative error in three-term recurrence:", err
            call exit(1)
        endif
    enddo

    if (verbose) then
        print *, "o Three-term recurrence for the Lanczos process was"
        print *, "  correct up to machine precision!"
    endif

    ! Check that the Lanczos vectors are orthogonal by forming the matrix
    ! Q = V^T * V; `Q` should be about equal to the identity matrix.
    Q = matmul(transpose(V), V)

    ! Subtract the identity matrix from `Q`.
    do i = 1, nq
        Q(i, i) = Q(i, i) - 1.0_dp
    enddo

    ! Compute the Frobenius norm of `Q`
    Q = matmul(transpose(Q), Q)

    err = 0.0_dp
    do i = 1, nq
       err = err + Q(i, i)
    enddo
    err = dsqrt(err) / nq

    if (err > 1.0e-14) then
        print *, "Lanczos vectors are not orthogonal!"
        print *, "|| V^t * V - I ||_F = ", err
        call exit(1)
    endif

    if (verbose) then
        print *, "o Lanczos vectors are orthogonal to machine precision!"
    endif

end program eigensolver_test_lanczos
