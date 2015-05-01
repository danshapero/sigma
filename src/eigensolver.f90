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
subroutine lanczos(A, T, Q)                                                !
!--------------------------------------------------------------------------!
! Perform n steps of the Lanczos process on the matrix A, putting the      !
! coefficients into the 3 x n array T and the vectors into the array Q.    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    real(dp), intent(out) :: T(:,:), Q(:,:)
    ! local variables
    integer :: i, k, n
    real(dp) :: w(A%nrow), alpha, beta

    n = size(T, 2)

    T = 0.0_dp
    Q = 0.0_dp
    w = 0.0_dp

    !-----------------------------------
    ! First step of the Lanczos process
    call init_seed()

    ! Make a random unit vector
    call random_number(Q(:,1))
    Q(:,1) = 2 * Q(:,1) - 1
    Q(:,1) = Q(:,1) / dsqrt(sum(Q(:,1) * Q(:,1)))

    ! Fill out the entries in T, Q
    call A%matvec(Q(:,1), w)
    alpha = sum(Q(:,1) * w)
    w = w - alpha * Q(:,1)
    beta = dsqrt(sum(w * w))
    Q(:,2) = w / beta
    T(2,1) = alpha
    T(3,1) = beta
    T(1,1) = beta


    !-----------------------------------------
    ! Subsequent steps of the Lanczos process
    do i = 2, n - 1
        call A%matvec(Q(:,i), w)
        alpha = sum(Q(:,i) * w)
        w = w - alpha * Q(:,i) - beta * Q(:,i - 1)

        ! Re-orthogonalize the current vector against previous Lanczos vectors
        ! to avoid floating point goofs.
        do k = 1, i - 2
            w = w - sum(Q(:,k) * w) * Q(:,k)
        enddo

        beta = dsqrt(sum(w * w))
        Q(:,i+1) = w / beta
        T(2,i) = alpha
        T(3,i) = beta
        T(1,i) = beta
    enddo

    !---------------------------------------
    ! Last iteration of the Lanczos process
    call A%matvec(Q(:,n) , w)
    T(2, n) = sum(Q(:,n) * w)

end subroutine lanczos



!--------------------------------------------------------------------------!
subroutine generalized_lanczos(A, B, T, Q)                                 !
!--------------------------------------------------------------------------!
! Same as above routine, only for generalized eigenvalue problem           !
!     A * x = lambda * B * x.                                              !
! Note: every step of the Lanczos process for the generalized eigenvalue   !
! problem necessitates a solution of the system B * x = y. It is assumed   !
! that `B` has a solver for it set.                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A, B
    real(dp), intent(out) :: T(:,:), Q(:,:)
    ! local variables
    integer :: i, k, n
    real(dp) :: w(A%nrow), v(A%nrow), alpha, beta
    real(dp), allocatable :: z(:,:)

    n = size(T, 2)
    allocate(z(A%nrow, 0:n))

    T = 0.0_dp
    Q = 0.0_dp
    w = 0.0_dp
    v = 0.0_dp
    z = 0.0_dp

    !-----------------------------------
    ! First step of the Lanczos process
    call init_seed()

    ! Make a random unit vector
    call random_number(Q(:,1))
    Q(:,1) = 2 * Q(:,1) - 1
    call B%matvec(Q(:,1), w)
    Q(:,1) = Q(:,1) / dsqrt(sum(w * Q(:,1)))
    call B%matvec(Q(:,1), z(:,1))

    alpha = 0.0_dp
    beta = 0.0_dp

    do i = 1, n - 1
        call A%matvec(Q(:,i), w)
        v = w - beta * z(:, i-1)
        alpha = sum(v * Q(:,i))
        v = v - alpha * z(:,i)

        call B%solve(w, v)

        beta = dsqrt(sum(w * v))
        Q(:, i+1) = w / beta
        z(:, i+1) = v / beta

        T(2, i) = alpha
        T(3, i) = beta
        T(1, i) = beta
    enddo

    call A%matvec(Q(:,n), v)
    v = v - beta * z(:,n)
    T(2, n) = sum(Q(:,i) * v)

end subroutine generalized_lanczos



!--------------------------------------------------------------------------!
subroutine eigensolve(A, lambda, V)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    real(dp), intent(out) :: lambda(:), V(:,:)
    ! local variables
    integer :: i, info, n
    real(dp), allocatable :: T(:,:), Q(:,:), work(:)

    n = size(lambda)
    allocate(T(3, n), Q(n, n), work(2 * n - 2))

    call lanczos(A, T, V)

    call dstev('V', n, T(2, 1:n), T(3, 1:n-1), Q, n, work, info)

    V = matmul(V, Q)

    do i = 1, n
        V(:,i) = V(1,i) / dabs(V(1, i)) * V(:,i)
    enddo

    lambda = T(2, 1:n)

end subroutine eigensolve



!--------------------------------------------------------------------------!
subroutine generalized_eigensolve(A, B, lambda, V)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A, B
    real(dp), intent(out) :: lambda(:), V(:,:)
    ! local variables
    integer :: i, info, n
    real(dp), allocatable :: T(:,:), Q(:,:), work(:)

    n = size(lambda)
    allocate(T(3, n), Q(n, n), work(2 * n - 2))

    call generalized_lanczos(A, B, T, V)
    call dstev('V', n, T(2, 1:n), T(3, 1:n-1), Q, n, work, info)

    V = matmul(V, Q)

    lambda = T(2, 1:n)

end subroutine generalized_eigensolve


end module eigensolver
