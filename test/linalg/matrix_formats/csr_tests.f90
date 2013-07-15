program csr_tests

    use csr
    use matrix

    implicit none

    ! sparse matrix
    class(sparse_matrix), allocatable :: A

    ! variables for building the sparse matrix
    integer, allocatable :: rows(:),cols(:)
    real(kind(1d0)), allocatable :: vals(:)

    ! variables for testing accessors/mutators
    real(kind(1d0)) :: As(3,3),Bs(3,3)

    ! variables for testing matvec
    real(kind(1d0)), allocatable :: x(:),y(:)

    ! permutation


    ! some other locals
    integer :: i,j,k
    real(kind(1d0)) :: q



!--------------------------------------------------------------------------!
! Test building the matrix                                                 !
!--------------------------------------------------------------------------!
    allocate(csr_matrix::A)
    call A%init(100,100,298)

    allocate(rows(298),cols(298))
    do i=1,99
        rows(i) = i
        cols(i) = i

        rows(i+100) = i
        cols(i+100) = i+1

        rows(i+199) = i+1
        cols(i+199) = i
    enddo

    rows(100) = 100
    cols(100) = 100

    call A%build(rows,cols)

    do i=1,99
        call A%set_value(i,i,2.d0)
        call A%set_value(i,i+1,-1.d0)
        call A%set_value(i+1,i,-1.d0)
    enddo
    call A%set_value(100,100,2.d0)

    print *, A%get_value(1,100),A%nnz
    call A%set_value(1,100,-1.d0)
    print *, A%get_value(1,100),A%nnz

end program csr_tests
