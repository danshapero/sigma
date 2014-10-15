!==========================================================================!
module sigma                                                               !
!==========================================================================!
!==== This is a super-module for all of SiGMA; use it to include the   ====!
!==== entire library.                                                  ====!
!==========================================================================!

use graphs
use linear_operators
use sparse_matrices

! Use the solver and preconditioner modules
use cg_solvers
use bicgstab_solvers
use jacobi_solvers
use ldu_solvers

!use wrapper

! Use other auxiliary modules
use eigensolver
use vectors




end module sigma
