program matrix_tests

use fempack

implicit none

    class(graph), pointer :: g
    type(sparse_matrix) :: A
    integer :: i,j,k,frmt,ordering
    real(dp) :: z
    logical :: correct
    integer, allocatable :: edges(:,:), nbrs(:), p(:)
    real(dp), allocatable :: x(:), y(:)
    character(len=3) :: orientation

    allocate(edges(2,31), x(7), y(7), nbrs(8))
    edges(:,1)  = [1, 2]
    edges(:,2)  = [1, 3]
    edges(:,3)  = [1, 4]
    edges(:,4)  = [1, 5]
    edges(:,5)  = [1, 6]
    edges(:,6)  = [1, 7]
    edges(:,7)  = [2, 3]
    edges(:,8)  = [3, 4]
    edges(:,9)  = [4, 5]
    edges(:,10) = [5, 6]
    edges(:,11) = [6, 7]
    edges(:,12) = [7, 2]
    do k=1,12
        edges(1,k+12) = edges(2,k)
        edges(2,k+12) = edges(1,k)
    enddo
    do k=1,7
        edges(:,k+24) = [k,k]
    enddo


    ! Loop through all formats
    do frmt=1,4
        do ordering=1,2
            ! Choose whether the graph is to be stored in row- or column-
            ! major ordering
            if (ordering==1) then
                orientation = 'row'
            else
                orientation = 'col'
            endif

            ! Choose the matrix format
            select case(frmt)
                case(1)
                    allocate(ll_graph::g)
                case(2)
                    allocate(coo_graph::g)
                case(3)
                    allocate(cs_graph::g)
                case(4)
                    allocate(ellpack_graph::g)
            end select
            call g%init(7,7,edges)

            ! Build the matrix
            call A%init(7,7,orientation,g)

            do i=1,7
                call A%g%neighbors(i,nbrs)

                do k=1,A%max_degree
                    j = nbrs(k)
                    if (j/=0 .and. j/=i) then
                        call A%set_value(i,j,-1.0_dp)
                        call A%add_value(i,i,1.0_dp)
                    endif
                enddo
            enddo

            ! Check that the matrix values were set right
            z = A%get_value(1,2)
            if (z/=-1.0_dp) then
                print *, 'On matrix test',frmt,ordering
                print *, 'A(1,2) should be = -1.0'
                print *, 'Value found: ',z
                call exit(1)
            endif

            z = A%get_value(1,1)
            if (z/=6.0_dp) then
                print *, 'On matrix test',frmt,ordering
                print *, 'A(1,1) should be = 6.0'
                print *, 'Value found: ',z
                call exit(1)
            endif

            ! Test matrix multiplication
            y = 2.0_dp
            x = 1.0_dp

            call A%matvec_add(x,y)
            if (minval(dabs(y))/=2.0_dp .or. maxval(dabs(y))/=2.0_dp) then
                print *, 'On matrix test',frmt,ordering
                print *, 'Graph Laplacian * constant vector should = 0.0'
                print *, 'Value found:',maxval(dabs(y))-2.0_dp
            endif



            ! Deallocate the graph and free the matrix
            deallocate(g)
            call A%destroy()
        enddo
    enddo

    deallocate(x,y,nbrs,edges)

end program matrix_tests
