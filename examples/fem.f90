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
    integer :: i,j,n,ne
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



!--------------------------------------------------------------------------!
!subroutine system2d(A,x,ele,kappa,d)                                       !
!--------------------------------------------------------------------------!
    ! input/output variables
!    class(block_sparse_matrix), intent(inout) :: A
!    real(dp), intent(in) :: x(:,:), kappa(:,:,:,:,:)
!    integer, intent(in) :: ele(:,:), d
    ! local variables
!    integer :: i,j,k,l,n,ne,elem(3)
!    real(dp) :: det,area,S(2,2),T(2,2),V(3,2),grad(3,2),kap(d,2,d,2), &
!        & AE(d,3,d,3)

!    ne = size(ele,2)
    ! call A%zero()

!    V(1,:) = [ 1.0_dp,  0.0_dp ]
!    V(2,:) = [ 0.0_dp,  1.0_dp ]
!    V(3,:) = [-1.0_dp, -1.0_dp ]

!    do n=1,ne
!        elem = ele(:,n)

        ! y -> x_3 + T*y maps reference triangle to physical triangle
!        T(:,1) = x(:,elem(1))-x(:,elem(3))
!        T(:,2) = x(:,elem(2))-x(:,elem(3))
!        det = T(1,1)*T(2,2)-T(1,2)*T(2,1)
!        area = 0.5*dabs(det)

        ! x -> S*(x-x_3) maps physical triangle to reference triangle
!        S(1,1) =  T(2,2)/det
!        S(1,2) = -T(1,2)/det
!        S(2,1) = -T(2,1)/det
!        S(2,2) =  T(1,1)/det

        ! Compute the gradient of all the shape functions
!        grad = matmul(V,S)

        ! Do something less blunt than averaging the diffusivity tensor
!        kap = sum( kappa(:,:,:,:,elem),5 )/3.0_dp

        ! Fill in the entries of the element stiffness matrix
!        AE = 0.0_dp
!        do l=1,d
!            do k=1,d
!                AE(k,:,l,:) = &
!                    &+area*matmul(grad,matmul(kap(k,:,l,:),transpose(grad)))
!            enddo
!        enddo

        ! Add the element stiffness matrix to the global matrix
!        do j=1,3
!            do i=1,3
!                call A%add_block(elem(i),elem(j),AE(:,i,:,j))
!            enddo
!        enddo

!    enddo

!end subroutine system2d



!--------------------------------------------------------------------------!
subroutine gradient(x,ele,u,grad)                                          !
!--------------------------------------------------------------------------!
    ! input/output variables
    real(dp), intent(in) :: x(:,:), u(:)
    integer, intent(in) :: ele(:,:)
    real(dp), intent(out) :: grad(:,:)
    ! local variables
    integer :: n,ne,elem(3)
    real(dp) :: det,S(2,2),T(2,2)

    ne = size(ele,2)

    do n=1,ne
        elem = ele(:,n)

        ! y -> x_3 + T*y maps reference triangle to physical triangle
        T(:,1) = x(:,elem(1))-x(:,elem(3))
        T(:,2) = x(:,elem(2))-x(:,elem(3))

        ! x -> S*(x-x_3) maps physical triangle to reference triangle
        det = T(1,1)*T(2,2)-T(1,2)*T(2,1)
        S(1,1) =  T(2,2)/det
        S(2,2) =  T(1,1)/det
        S(2,1) = -T(2,1)/det
        S(1,2) = -T(1,2)/det

        ! Compute the directional derivatives of u along the edges of the
        ! physical triangle
        grad(1,n) = u(elem(1))-u(elem(3))
        grad(2,n) = u(elem(2))-u(elem(3))

        ! Multiply by S to find the derivatives in the x_1, x_2 directions
        grad(:,n) = matmul(grad(:,n),S)
    enddo

end subroutine gradient





end module fem
