!--------------------------------------------------------------------------!
program matrix_test_ptap                                                   !
!--------------------------------------------------------------------------!
! This program tests explicitly forming the product P^t * A * P of the     !
! sparse matrices A, P rather than lazily forming the operator product.    !
!--------------------------------------------------------------------------!

use sigma

implicit none

    ! graphs used as the matrix substrates
    type(ll_graph) :: g, h

    ! sparse and dense matrices
    class(sparse_matrix_interface), pointer :: A, B, P
    real(dp), allocatable :: AD(:,:), BD(:,:), PD(:,:)

    ! indices
    integer :: i, j, k, l, n, nn, nnc

    ! variables for getting graph neighbors
    integer :: d
    integer, allocatable :: nodes(:)

    ! error in computing sparse matrix product
    real(dp) :: misfit

    ! command-line argument parsing
    character(len=16) :: arg
    logical :: verbose

    ! random numbers
    real(dp) :: z, probability



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

    call init_seed()


    !----------------------------------------------------------------------!
    ! Make a random graph and a coarsened approximation to it              !
    !----------------------------------------------------------------------!

    nn = 256
    call g%init(nn)

    probability = 0.5 * log(1.0_dp * nn) / log(2.0_dp) / nn
    do i = 1, nn
        do j = i + 1, nn
            call random_number(z)

            if (z < probability) then
                call g%add_edge(i, j)
                call g%add_edge(j, i)
            endif
        enddo

        call g%add_edge(i, i)
    enddo

    if (verbose) then
        print *, "o Done generating random graph."
        print *, "    Number of vertices:", nn
        print *, "    Number of edges:   ", g%get_num_edges()
        print *, "    Max degree:        ", g%get_max_degree()
    endif


    nnc = 128
    call h%init(nn, nnc)

    do i = 1, nn
        do j = 1, nnc
            call random_number(z)

            if (z < probability / 2) then
                call h%add_edge(i, j)
            endif
        enddo
    enddo

    if (verbose) then
        print *, "o Done generating graph of prolongation matrix."
        print *, "    Number of vertices:", nn, nnc
        print *, "    Number of edges:   ", h%get_num_edges()
        print *, "    Max degree:        ", h%get_max_degree()
    endif


    d = max(g%get_max_degree(), h%get_max_degree())
    allocate(nodes(d))



    !----------------------------------------------------------------------!
    ! Make sparse matrices on those graphs                                 !
    !----------------------------------------------------------------------!

    allocate(csc_matrix :: A)
    allocate(csr_matrix :: P)

    call A%init(nn, nn, g)
    call P%init(nn, nnc, h)

    call A%zero()
    call P%zero()

    do j = 1, nn
        d = g%get_degree(j)
        call g%get_neighbors(nodes, j)

        do k = 1, d
            i = nodes(k)

            call A%add(j, j, +1.0_dp)
            call A%add(i, j, -1.0_dp)
        enddo
    enddo


    do i = 1, nn
        d = h%get_degree(i)
        call h%get_neighbors(nodes, i)

        do k = 1, d
            j = nodes(k)
            call P%set(i, j, 1.0_dp / d)
        enddo
    enddo



    !----------------------------------------------------------------------!
    ! Compute congruent product B = P^t * A * P                            !
    !----------------------------------------------------------------------!

    allocate(csc_matrix :: B)
    call PtAP(B, A, P)


    if (verbose) then
        print *, "o Done computing product P^t * A * P."
        print *, "    Number of vertices:", B%nrow, B%ncol
        print *, "    Number of edges:   ", B%get_nnz()
        print *, "    Max degree:        ", B%get_max_column_degree()
    endif



    !----------------------------------------------------------------------!
    ! Check that computing the congruent product worked                    !
    !----------------------------------------------------------------------!

    allocate(AD(nn, nn), PD(nn, nnc), BD(nnc, nnc))

    AD = 0.0_dp
    PD = 0.0_dp
    BD = 0.0_dp

    call A%to_dense_matrix(AD)
    call P%to_dense_matrix(PD)
    call B%to_dense_matrix(BD)


    misfit = maxval(dabs( BD - matmul(transpose(PD), matmul(AD, PD)) ))


    if (misfit > 1.0e-15) then
        print *, "Computing congruent product B = P^t * A * P failed."
        print *, "Misfit:", misfit
        call exit(1)
    endif


    call A%destroy()
    call P%destroy()
    call B%destroy()

    deallocate(A)
    deallocate(P)
    deallocate(B)


end program matrix_test_ptap

