program self_avoiding_walk

use sigma
use regular_graphs
use random_graphs

implicit none

    ! variables for parsing command line arguments
    character(len=32) :: arg, graph_name, int_param1, int_param2, &
        & real_param, iter_name
    logical :: verbose

    ! variables for generating a graph which will be the "medium" of the
    ! self-avoiding walk
    class(graph), pointer :: g
    integer :: n, k
    integer, allocatable :: s(:), neighbors(:), histogram(:)
    logical, allocatable :: unvisited(:)
    real(dp) :: beta, p

    ! other variables
    integer :: i,j,d,nn,iter,trial,num
    real(dp) :: z
    type(circular_array) :: q

    graph_name = "torus                           "
    int_param1 = "32                              "
    int_param2 = "4                               "
    real_param = "0.25                            "
    iter_name =  "10000                           "
    verbose = .false.

    !----------------------------------------------------------------------!
    ! Get command-line arguments                                           !
    !----------------------------------------------------------------------!
    do i=1,iargc()
        call getarg(i,arg)
        select case(trim(arg))
            case("--graph","-g")
                call getarg(i+1,graph_name)
            case("--n","-n")
                call getarg(i+1,int_param1)
            case("--k","-k")
                call getarg(i+1,int_param2)
            case("--p","-p")
                call getarg(i+1,real_param)
            case("--iter","-i")
                call getarg(i+1,iter_name)
            case("-v","--verbose","-V")
                verbose = .true.
        end select
    enddo

    read(int_param1,*) n
    read(int_param2,*) k
    read(real_param,*) p
    read(iter_name,*) iter


    !----------------------------------------------------------------------!
    ! Create the "medium", a graph                                         !
    !----------------------------------------------------------------------!
    select case(trim(graph_name))
        case("torus")
            allocate(cs_graph::g)
            if (verbose) then
                write(*,10) n,k
10              format('Constructing a 2-torus of dimensions ',i4,' x ',i4)
            endif
            call torus(g,n,k)


        case("petersen")
            allocate(cs_graph::g)
            if (verbose) then
                write(*,20) n,k
20              format('Constructing a Petersen graph with ',i4,',',i4)
            endif
            call petersen(g,n,k)


        case("snark","flower-snark","flower_snark","flowersnark")
            allocate(cs_graph::g)
            if (verbose) then
                write(*,30) 4*n
30              format('Constructing a flower snark on ',i4,' vertices.')
            endif
            call flower_snark(g,n)


        case("hypercube")
            allocate(cs_graph::g)
            k = min(n,10)
            if (verbose) then
                write(*,40) 2**k
40              format('Constructing hypercube graph on ',i4,' vertices.')
            endif
            call hypercube(g,k)


        case("erdos-renyi","erdos_renyi","erdosrenyi","er")
            allocate(ll_graph::g)
            call erdos_renyi(g,n,(1.0_dp*k)/n)
        case("watts-strogatz","watts_strogatz","wattsstrogatz","ws", &
                & "small-world","small_world","smallworld")
            allocate(ll_graph::g)
            call watts_strogatz(g,n,k,p)
        case("barabasi-albert","barabasi_albert","barabasialbert","ba", &
                & "scale-free","scale_free","scalefree")
            allocate(ll_graph::g)
            call barabasi_albert(g,n,k)
    end select

    call g%compress()
    nn = g%n
    allocate(neighbors(g%max_degree), unvisited(g%max_degree), s(nn))


    !----------------------------------------------------------------------!
    ! Execute a bunch of self-avoiding walks                               !
    !----------------------------------------------------------------------!
    call init_seed()
    allocate(histogram(0:nn))
    histogram = 0

    do trial=1,iter
        ! Choose a random starting vertex
        call random_number(z)
        i = min(int(z*nn)+1,nn)

        ! Initialize the list indicating which vertices have been visited
        ! to zero
        s = 0

        ! Initialize a queue of visited vertices
        call q%init(capacity=g%max_degree,min_capacity=2)

        ! So long as there are unvisited neighbor vertices, keep walking
        unvisited = .true.
        do while (any(unvisited))
            ! Put the next vertex in the queue and indicate that it's been
            ! visited
            call q%enqueue(i)
            s(i) = q%length

            ! Assume that all the neighbors of the vertex i have been visited
            unvisited = .false.

            ! Check to see which neighbors of i have not been visited
            d = g%degree(i)
            call g%get_neighbors(neighbors,i)

            do k=1,d
                j = neighbors(k)
                unvisited(k) = ( s(j)==0 .and. j/=i )
            enddo

            ! Pick one of the unvisited neighbors at random, then make it
            ! the next vertex
            call random_number(z)
            num = min(int(count(unvisited)*z)+1,d)

            do k=1,d
                j = neighbors(k)

                if (unvisited(k) .and. num>0) then
                    num = num-1
                    i = j
                endif
            enddo
        enddo

        num = q%length
        histogram(num) = histogram(num)+1

        ! Clear the queue
        call q%free()
    enddo

    do i=0,nn
        print *, histogram(i)
    enddo


end program self_avoiding_walk
