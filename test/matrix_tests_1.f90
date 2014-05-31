!--------------------------------------------------------------------------!
program matrix_tests_1                                                     !
!--------------------------------------------------------------------------!
!     This program performs tests of basic matrix operations, namely:      !
! initialization, getting and setting entries, matrix multiplication,      !
! storage compression and permutation.                                     !
!--------------------------------------------------------------------------!

use sigma

implicit none

    class(graph), pointer :: g
    type(sparse_matrix) :: A
    integer :: i,j,k,d,frmt,ordering
    real(dp) :: z
    logical :: correct
    integer, allocatable :: edges(:,:), neighbors(:), p(:)
    real(dp), allocatable :: x(:), y(:), vals(:)
    character(len=3) :: orientation
    real(dp), allocatable :: B(:,:)

    allocate(edges(2,31), x(7), y(7), neighbors(8), vals(8), B(7,7))

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
            call choose_graph_type(g,frmt)
            call g%init(7,7,[7,4,4,4,4,4,4])

            do i=1,6
                call g%add_edge(i,i+1)
                call g%add_edge(i+1,i)
                call g%add_edge(1,i+1)
                call g%add_edge(i+1,1)
                call g%add_edge(i,i)
            enddo
            call g%add_edge(7,2)
            call g%add_edge(2,7)
            call g%add_edge(7,7)

            ! Build the matrix
            call A%init(7,7,orientation,g)

            do i=1,7
                call A%g%get_neighbors(neighbors,i)

                do k=1,A%max_degree
                    j = neighbors(k)
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

            ! Check that getting a matrix row/column works
            do i=1,7
                d = A%g%degree(i)
                call A%get_row(neighbors,vals,i)

                do k=1,d
                    j = neighbors(k)
                    z = A%get_value(i,j)

                    if (vals(k)/=z) then
                        print *, 'On matrix test',frmt,ordering
                        print *, 'Getting row',i,'failed;'
                        print *, 'Should have entry',i,j,' = ',z
                        print *, 'Value found:',vals(k)
                        print *, 'Terminating.'
                        call exit(1)
                    endif
                enddo
            enddo

            do j=1,7
                d = A%g%degree(j)
                call A%get_column(neighbors,vals,j)

                do k=1,d
                    i = neighbors(k)
                    z = A%get_value(i,j)

                    if (vals(k)/=z) then
                        print *, 'On matrix test',frmt,ordering
                        print *, 'Getting column',j,'failed;'
                        print *, 'Should have entry',i,j,' = ',z
                        print *, 'Value found:',vals(k)
                        print *, 'Terminating.'
                        call exit(1)
                    endif
                enddo
            enddo

            ! Test matrix multiplication
            y = 2.0_dp
            x = 1.0_dp

            call A%matvec_add(x,y)
            if (minval(dabs(y))/=2.0_dp .or. maxval(dabs(y))/=2.0_dp) then
                print *, 'On matrix test',frmt,ordering
                print *, 'Graph Laplacian * constant vector should = 0.0'
                print *, 'Value found:',maxval(dabs(y))-2.0_dp
                call exit(1)
            endif

            ! Compress the matrix and test to see if matrix multiplication
            ! still works
            call A%compress()

            call A%matvec_add(x,y)

            ! Set all the matrix entries so it's easiest to tell if we've
            ! permuted everything right, and set up an equivalent dense
            ! matrix to check against
            B = 0.0_dp
            do i=1,7
                call A%g%get_neighbors(neighbors,i)
                do k=1,A%max_degree
                    j = neighbors(k)
                    if (j/=0) then
                        call A%set_value(i,j,1.0_dp*(7*(i-1)+j))
                        B(i,j) = 1.0_dp*(7*(i-1)+j)
                    endif
                enddo
            enddo

            ! Make a permutation (7 1 2 3 4 5 6)
            allocate(p(7))
            p(1) = 7
            do i=2,7
                p(i) = i-1
            enddo

            ! Check that right-permutation works correctly
            B(:,p) = B(:,:)
            call A%right_permute(p)
            correct = .true.
            do j=1,7
                do i=1,7
                    correct = correct .and. (A%get_value(i,j)==B(i,j))
                enddo
            enddo
            if (.not.correct) then
                print *, 'On matrix test',frmt,ordering
                print *, 'Permuting the columns of the matrix failed'
                call exit(1)
            endif

            ! Check that left-permutation works properly
            B(p,:) = B(:,:)
            call A%left_permute(p)
            correct = .true.
            do j=1,7
                do i=1,7
                    correct = correct .and. (A%get_value(i,j)==B(i,j))
                enddo
            enddo
            if (.not.correct) then
                print *, 'On matrix test',frmt,ordering
                print *, 'Permuting the rows of the matrix failed'
                call exit(1)
            endif

            ! Deallocate the graph and destroy the matrix
            call g%destroy()
            deallocate(g,p)
            call A%destroy()
        enddo
    enddo

    deallocate(x,y,neighbors,edges)

end program matrix_tests_1
