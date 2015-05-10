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
subroutine lanczos(A, alpha, beta, Q)                                      !
!--------------------------------------------------------------------------!
! Perform n steps of the Lanczos process on the matrix A, putting the      !
! coefficients into the 3 x n array T and the vectors into the array Q.    !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    real(dp), intent(out) :: alpha(:), beta(:), Q(:,:)
    ! local variables
    integer :: i, k, n
    real(dp) :: w(A%nrow)

    n = size(alpha)

    alpha = 0.0_dp
    beta  = 0.0_dp
    Q     = 0.0_dp
    w     = 0.0_dp

    !-----------------------------------
    ! First step of the Lanczos process
    call init_seed()

    ! Make a random unit vector
    call random_number(Q(:,1))
    Q(:,1) = 2 * Q(:,1) - 1
    Q(:,1) = Q(:,1) / dsqrt(dot_product(Q(:,1), Q(:,1)))

    call A%matvec(Q(:,1), w)
    alpha(1) = dot_product(w, Q(:, 1))
    w = w - alpha(1) * Q(:, 1)
    beta(2) = dsqrt(dot_product(w, w))
    Q(:,2) = w / beta(2)

    !-----------------------------------------
    ! Subsequent steps of the Lanczos process
    do i = 2, n - 1
        call A%matvec(Q(:,i), w)
        w = w - beta(i) * Q(:,i-1)
        alpha(i) = dot_product(w, Q(:,i))
        w = w - alpha(i) * Q(:,i)

        ! Re-orthogonalize the current vector against previous Lanczos vectors
        ! to avoid floating point goofs.
        do k = 1, i - 2
            w = w - sum(Q(:,k) * w) * Q(:,k)
        enddo

        beta(i + 1) = dsqrt(dot_product(w, w))
        Q(:,i+1) = w / beta(i + 1)
    enddo

    !---------------------------------------
    ! Last iteration of the Lanczos process
    call A%matvec(Q(:,n), w)
    w = w - beta(n) * Q(:, n - 1)
    alpha(n) = dot_product(Q(:,n), w)

end subroutine lanczos



!--------------------------------------------------------------------------!
subroutine generalized_lanczos(A, B, alpha, beta, Q)                       !
!--------------------------------------------------------------------------!
! Same as above routine, only for generalized eigenvalue problem           !
!     A * x = lambda * B * x.                                              !
! Note: every step of the Lanczos process for the generalized eigenvalue   !
! problem necessitates a solution of the system B * x = y. It is assumed   !
! that `B` has a solver for it set.                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A, B
    real(dp), intent(out) :: alpha(:), beta(:), Q(:,:)
    ! local variables
    integer :: i, k, n
    real(dp) :: w(A%nrow), v(A%nrow), r
    real(dp), allocatable :: z(:,:)

    n = size(alpha)
    allocate(z(A%nrow, 0:n))

    alpha = 0.0_dp
    beta = 0.0_dp
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
    Q(:,1) = Q(:,1) / dsqrt(dot_product(w, Q(:,1)))
    call B%matvec(Q(:,1), z(:,1))

    do i = 1, n - 1
        call A%matvec(Q(:,i), v)
        v = v - beta(i) * z(:, i-1)
        alpha(i) = dot_product(v, Q(:,i))
        v = v - alpha(i) * z(:,i)

        call B%solve(w, v)

        beta(i + 1) = dsqrt(dot_product(w, v))
        Q(:, i+1) = w / beta(i + 1)
        z(:, i+1) = v / beta(i + 1)
    enddo

    call A%matvec(Q(:,n), v)
    v = v - beta(n) * z(:,n-1)
    alpha(n) = dot_product(Q(:,i), v)

end subroutine generalized_lanczos



!--------------------------------------------------------------------------!
subroutine eigensolve(A, lambda, V)                                        !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A
    real(dp), intent(out) :: lambda(:), V(:,:)
    ! local variables
    integer :: i, info, n
    real(dp), allocatable :: alpha(:), beta(:), Q(:,:), work(:)

    n = size(lambda)
    allocate(alpha(n), beta(n), Q(n, n), work(2 * n - 2))

    call lanczos(A, alpha, beta, V)

    call dstev('V', n, alpha, beta(2:n), Q, n, work, info)

    V = matmul(V, Q)

    do i = 1, n
        V(:,i) = V(1,i) / dabs(V(1, i)) * V(:,i)
    enddo

    lambda = alpha

end subroutine eigensolve



!--------------------------------------------------------------------------!
subroutine generalized_eigensolve(A, B, lambda, V)                         !
!--------------------------------------------------------------------------!
    ! input/output variables
    class(sparse_matrix_interface), intent(in) :: A, B
    real(dp), intent(out) :: lambda(:), V(:,:)
    ! local variables
    integer :: i, info, n
    real(dp), allocatable :: alpha(:), beta(:), Q(:,:), work(:)

    n = size(lambda)
    allocate(alpha(n), beta(n), Q(n, n), work(2 * n - 2))

    call generalized_lanczos(A, B, alpha, beta, V)
    call dstev('V', n, alpha, beta(2:n), Q, n, work, info)

    V = matmul(V, Q)

    lambda = alpha

end subroutine generalized_eigensolve


end module eigensolver
