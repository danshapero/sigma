module fem

use fempack

implicit none

contains


!--------------------------------------------------------------------------!
subroutine laplacian2d(A,x,ele)                                            !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: A
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: ele(:,:)
    ! local variables
    integer :: i,j,k,n,ne
    real(dp) :: AE(3,3), V(3,2), det, area

    ne = size(ele,2)

    if (size(x,1)/=2) then
        print *, 'Wrong dimension for input array x'
        call exit(1)
    endif

    do n=1,ne
        AE = 0.0_dp

        do i=1,3
            j = ele(mod(i,3)+1,n)
            k = ele(mod(i+1,3)+1,n)
            V(i,1) = x(2,j)-x(2,k)
            V(i,2) = x(1,k)-x(1,j)
        enddo

        det = V(1,1)*V(2,2)-V(1,2)*V(2,1)
        area = dabs(det)/2.0_dp

        AE = 0.25/area*matmul(V,transpose(V))

        do j=1,3
            do i=1,3
                call A%add_value(ele(i,n),ele(j,n),AE(i,j))
            enddo
        enddo

    enddo

end subroutine laplacian2d



!--------------------------------------------------------------------------!
subroutine mass2d(B,x,ele)                                                 !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(inout) :: B
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: ele(:,:)
    ! local variables
    integer :: i,j,k,n,ne
    real(dp) :: BE(3,3), area

    ne = size(ele,2)

    do n=1,ne
        BE = 0.0_dp
        do j=1,2
            do i=1,2
                BE(i,j) = x(i,ele(j,n))-x(i,ele(3,n))
            enddo
        enddo
        area = 0.5*dabs(BE(1,1)*BE(2,2)-BE(1,2)*BE(2,1))

        BE = area/12.0_dp
        do i=1,3
            BE(i,i) = area/6.0_dp
        enddo

        do j=1,3
            do i=1,3
                call B%add_value(ele(i,n),ele(j,n),BE(i,j))
            enddo
        enddo
    enddo

end subroutine mass2d




end module fem
