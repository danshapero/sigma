program ising

use sigma
use regular_graphs
use random_graphs

implicit none

    ! variables for parsing command line arguments
    character(len=32) :: arg, graph_name, int_param1, int_param2, &
        & real_param, beta_name
    logical :: verbose

    ! variables for generating a graph which will be the "medium" of the
    ! Ising model
    class(graph), pointer :: g
    integer :: n, k
    integer, allocatable :: s(:), neighbors(:)
    real(dp) :: beta, p, dE

    ! variables for making a graph coloring
    !integer, allocatable :: color(:)

    ! other variables
    integer :: i, j, d, nn, iter
    real(dp) :: z, flip_probability
    real(dp), allocatable :: mag(:)
    real(dp) :: magnetization

    graph_name = "torus                           "
    int_param1 = "32                              "
    int_param2 = "4                               "
    real_param = "0.25                            "
    beta_name =  "1.0                             "
    verbose = .false.

    !----------------------------------------------------------------------!
    ! Get command-line arguments                                           !
    !----------------------------------------------------------------------!
    do i=1,iargc()
        call getarg(i,arg)
        select case(trim(arg))
            case("--graph", "-g")
                call getarg(i + 1, graph_name)
            case("--n", "-n")
                call getarg(i + 1, int_param1)
            case("--k", "-k")
                call getarg(i + 1, int_param2)
            case("--p", "-p")
                call getarg(i + 1, real_param)
            case("--beta", "-b", "inverseTemp", "inversetemp")
                call getarg(i + 1, beta_name)
            case("-v", "--verbose", "-V")
                verbose = .true.
        end select
    enddo

    read(int_param1,*) n
    read(int_param2,*) k
    read(real_param,*) p
    read(beta_name,*)  beta


    !----------------------------------------------------------------------!
    ! Create the "medium", a graph                                         !
    !----------------------------------------------------------------------!
    allocate(ll_graph::g)
    select case(trim(graph_name))
        case("torus")
            if (verbose) then
                write(*,10) n, k
10              format('Constructing a 2-torus of dimensions ',i4,' x ',i4)
            endif
            call torus(g, n, k)


        case("petersen")
            if (verbose) then
                write(*,20) n, k
20              format('Constructing a Petersen graph with ',i4,',',i4)
            endif
            call petersen(g, n, k)


        case("snark","flower-snark","flower_snark","flowersnark")
            if (verbose) then
                write(*,30) 4 * n
30              format('Constructing a flower snark on ',i4,' vertices.')
            endif
            call flower_snark(g, n)


        case("hypercube")
            k = min(n, 10)
            if (verbose) then
                write(*,40) 2**k
40              format('Constructing hypercube graph on ',i4,' vertices.')
            endif
            call hypercube(g, k)


        case("erdos-renyi","erdos_renyi","erdosrenyi","er")
            call erdos_renyi(g, n, (1.0_dp*k) / n)

        case("watts-strogatz","watts_strogatz","wattsstrogatz","ws", &
                & "small-world","small_world","smallworld")
            call watts_strogatz(g, n, k, p)

        case("barabasi-albert","barabasi_albert","barabasialbert","ba", &
                & "scale-free","scale_free","scalefree")
            call barabasi_albert(g, n, k)
    end select

    nn = g%n
    call convert_graph_type(g, 'compressed sparse')

    d = g%max_degree()
    allocate(neighbors(d))


    !----------------------------------------------------------------------!
    ! Find a coloring of the graph                                         !
    !----------------------------------------------------------------------!
!    allocate(color(nn))
!    call greedy_coloring(g,color)


    !----------------------------------------------------------------------!
    ! Make an initial, random configuration                                !
    !----------------------------------------------------------------------!
    allocate(s(nn))
    call init_seed()
    do i=1,nn
        !call random_number(z)
        !s(i) = 2*int(z+0.5)-1
        s(i) = 1
    enddo


    !----------------------------------------------------------------------!
    ! Execute several steps of the Metropolis algorithm                    !
    !----------------------------------------------------------------------!
    allocate(mag(nn))
    do iter = 1, 100 * nn
        do i = 1, nn
            call g%get_neighbors(neighbors, i)
            d = g%degree(i)
            do k = 1, d
                j = neighbors(k)

                dE = dE + s(j)
            enddo
            dE = dE * s(i)

            flip_probability = min(1.0_dp, exp(-beta * dE))
            call random_number(z)
            if (z <= flip_probability) s(i) = -s(i)
        enddo

        mag( mod(iter - 1, nn) + 1 ) = (1.0_dp*sum(s)) / nn

        if (mod(iter, nn)==0) then
            magnetization = (1.0_dp*sum(mag)) / nn
            print *, iter, magnetization
        endif
    enddo


    call g%destroy()


end program ising
