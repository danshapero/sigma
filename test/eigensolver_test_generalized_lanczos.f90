
!--------------------------------------------------------------------------!
program eigensolver_generalized_lanczos_test                               !
!--------------------------------------------------------------------------!
! Test the correctness of the Lanczos process for a generalized eigen-     !
! problem A*x = Lambda*B*x.                                                !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graph used as the matrix substrate
    type(ll_graph) :: g

    ! sparse matrices
    type(sparse_matrix) :: A, B

    ! element stiffness/mass matrices
    real(dp) :: AE(3, 3), BE(3, 3), area
    integer :: elem(3)

    ! ints
    integer :: i, j, k, l, ny, nx, nn, nq

    ! Lanczos vectors and tridiagonal matrix
    real(dp), allocatable :: U(:,:), V(:,:), T(:,:), Q(:,:)

    ! vectors
    real(dp), allocatable :: w(:), z(:)

    ! approximation error
    real(dp) :: err

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
    ! Make a graph representing a regular grid                             !
    !----------------------------------------------------------------------!

    nx = 48
    ny = 32
    nn = ny * nx

    call g%init(nn)

    do i = 1, ny
        do j = 1, nx
           k = indx(i, j)

           call g%add_edge(k, k)

           l = indx(mod(i, ny) + 1, j)
           call g%add_edge(k, l)
           call g%add_edge(l, k)

           l = indx(i, mod(j, nx) + 1)
           call g%add_edge(k, l)
           call g%add_edge(l, k)

           l = indx(mod(i, ny) + 1, mod(j, nx) + 1)
           call g%add_edge(k, l)
           call g%add_edge(l, k)
        enddo
    enddo

    if (verbose) then
        print *, "o Done creating graph of regular grid!"
        print *, "  nx: ", nx
        print *, "  ny: ", ny
        print *, "  number of neighbors:", g%get_max_degree()
    endif


    !----------------------------------------------------------------------!
    ! Create the stiffness & mass matrices on the grid                     !
    !----------------------------------------------------------------------!
    call A%set_matrix_type("csr")
    call B%set_matrix_type("csr")
    call A%set_dimensions(nn, nn)
    call B%set_dimensions(nn, nn)
    call A%copy_graph(g)
    call B%copy_graph(g)
    call A%zero()
    call B%zero()

    area = 0.5d0
    BE = area / 12.0_dp
    do k = 1, 3
        BE(k, k) = area / 6.0_dp
    enddo

    !    x
    !   /|
    !  / |
    ! x--x
    AE(:, 1) = [+area,  -area, 0.0_dp]
    AE(:, 2) = [-area,  2*area, -area]
    AE(:, 3) = [0.0_dp, -area,  +area]

    do i = 1, ny
        do j = 1, nx
            elem(1) = indx(i, j)
            elem(2) = indx(i, mod(j, nx) + 1)
            elem(3) = indx(mod(i, ny) + 1, mod(j, nx) + 1)

            call A%add(elem, elem, AE)
            call B%add(elem, elem, BE)

            elem(2) = indx(mod(i, ny) + 1, j)

            call A%add(elem, elem, AE)
            call B%add(elem, elem, BE)
        enddo
    enddo

    if (verbose) then
        print *, "o Done creating finite element matrices."
    endif


    !----------------------------------------------------------------------!
    ! Perform a few steps of the generalized Lanczos process               !
    !----------------------------------------------------------------------!
    nq = max(nx, ny)
    allocate(T(3, nq), U(nn, nq), V(nn, nq), Q(nq, nq), w(nn), z(nn))

    ! `B` must have a solver; it plays the role in the generalized Lanczos
    ! algorithm that the preconditioner does in PCG.
    ! Fortunately, finite-element mass matrices are usually well-conditioned
    ! if the underlying triangulation is nice enough.
    call B%set_solver(cg(1.0d-15))

    call generalized_lanczos(A, B, T, V)

    if (verbose) then
        print *, "o Done executing Lanczos algorithm."
    endif

    do i = 1, nq
        call B%matvec(V(:, i), U(:, i))
    enddo

    do i = 2, nq - 1
        call A%matvec(V(:, i), w)
        z = T(2, i) * U(:, i) + T(1, i-1) * U(:, i-1) + T(3, i) * U(:, i+1)

        err = dsqrt(sum((w-z) * (w-z)) / sum(w * w))

        if (err > 1.0e-14) then
            print *, "Computing Lanczos vector failed!"
            print *, "Relative error in three-term recurrence:", err
        endif
    enddo

    if (verbose) then
        print *, "o Three-term recurrence for the Lanczos process was"
        print *, "  correct up to machine precision!"
    endif

    ! Check that the Lanczos vectors are B-orthogonal by forming the matrix
    ! Q = V^T * B * V; `Q` should be about equal to the identity matrix.
    Q = matmul(transpose(V), U)

    ! Subtract the identity matrix from `Q`
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
        print *, "Lanczos vectors are not B-orthogonal!"
        print *, "|| V^t * B * V - I ||_F = ", err
    endif

    if (verbose) then
        print *, "o Lanczos vectors are B-orthogonal to machine precision!"
    endif

contains

    elemental function indx(i, j) result(k)
        integer, intent(in) :: i, j
        integer :: k

        k = ny * (j - 1) + i

    end function indx

    function coordinate(k) result(ij)
        integer, intent(in) :: k
        integer :: ij(2)

        ij(1) = mod(k - 1, ny) + 1
        ij(2) = (k - ij(1)) / ny + 1

    end function coordinate

end program eigensolver_generalized_lanczos_test
