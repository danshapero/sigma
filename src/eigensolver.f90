!==========================================================================!
!==========================================================================!
module eigensolver                                                         !
!==========================================================================!
!==========================================================================!
!==== This module contains subroutines for approximating eigenvalues   ====!
!==== and eigenvectors of linear operators through either the Lanczos  ====!
!==== or Arnoldi processes. Partial or full re-orthogonalization is    ====!
!==== used in order to stave off the breakdown of these algorithms.    ====!
!==========================================================================!
!==========================================================================!


use util, only: init_seed
use types, only: dp
use sparse_matrices

implicit none



contains



!--------------------------------------------------------------------------!
subroutine lanczos(A,T,Q,n)                                                !
!--------------------------------------------------------------------------!
! Perform n steps of the Lanczos process on the matrix A, putting the      !
! coefficients into the 3 x n array T and the vectors into the array Q.    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(out) :: T(3,n), Q(A%nrow,n)
    integer, intent(in) :: n
    ! local variables
    integer :: i,k
    real(dp) :: w(A%nrow), alpha, beta

    T = 0.0_dp
    Q = 0.0_dp
    w = 0.0_dp

    !-----------------------------------
    ! First step of the Lanczos process
    call init_seed()

    ! Make a random unit vector
    call random_number(Q(:,1))
    Q(:,1) = 2*Q(:,1)-1
    Q(:,1) = Q(:,1)/dsqrt(sum(Q(:,1)*Q(:,1)))

    ! Fill out the entries in T, Q
    call A%matvec(Q(:,1),w)
    alpha = sum(Q(:,1)*w)
    w = w-alpha*Q(:,1)
    beta = dsqrt(sum(w*w))
    Q(:,2) = w/beta
    T(2,1) = alpha
    T(3,1) = beta
    T(1,1) = beta
    

    !-----------------------------------------
    ! Subsequent steps of the Lanczos process
    do i=2,n-1
        call A%matvec(Q(:,i),w)
        alpha = sum(Q(:,i)*w)
        w = w-alpha*Q(:,i)-beta*Q(:,i-1)
        do k=1,i-2
            w = w-sum(Q(:,k)*w)*Q(:,k)
        enddo

        beta = dsqrt(sum(w*w))
        Q(:,i+1) = w/beta
        T(2,i) = alpha
        T(3,i) = beta
        T(1,i) = beta
    enddo

    !---------------------------------------
    ! Last iteration of the Lanczos process
    call A%matvec(Q(:,n),w)
    T(2,n) = sum(Q(:,n)*w)

end subroutine lanczos



!--------------------------------------------------------------------------!
subroutine eigensolve(A,lambda,V,n)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix), intent(in) :: A
    real(dp), intent(out) :: lambda(n), V(A%nrow,n)
    integer, intent(in) :: n
    ! local variables
    integer :: i,info
    real(dp) :: T(3,n), Q(n,n), work(2*n-2)

    call lanczos(A,T,V,n)

    call dstev('V',n,T(2,1:n),T(3,1:n-1),Q,n,work,info)

    V = matmul(V,Q)

    do i=1,n
        V(:,i) = V(1,i)/dabs(V(1,i))*V(:,i)
    enddo

    lambda = T(2,1:n)

end subroutine eigensolve


end module eigensolver
