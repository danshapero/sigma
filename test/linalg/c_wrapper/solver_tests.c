#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <matrix.h>
#include <solver.h>

int main(int argc, char *argv) {

    // Initialize the matrix
    sparse_matrix_c *A;
    A = (sparse_matrix_c *)malloc( sizeof(sparse_matrix_c) );
    get_sparse_matrix(A,0);

    int n = 99;
    init(A,n,n,3*n);

    // Build the non-zero structure of the matrix
    int rows[3*n], cols[3*n], i;
    /*for (i=0; i<n; i++) {
        rows[i] = i;
        rows[i+n] = i;
        rows[i+2*n] = i;
        cols[i] = i;
        cols[i+n] = i-1;
        cols[i+2*n] = i+1;
    }
    cols[n] = n-1;
    cols[2*n-1] = 0;

    build(A,rows,cols,3*n);

    // Fill in the matrix entries
    for (i=0; i<n-1; i++) {
        set_value(A,i,i,2.0);
        set_value(A,i,i+1,-1.0);
        set_value(A,i+1,i,-1.0);
    }
    set_value(A,n-1,n-1,2.0);
    set_value(A,0,n-1,0.0);
    set_value(A,n-1,0,0.0);


    // Set up the solver and preconditioner
    iterative_solver_c *solver;
    solver = (iterative_solver_c *)malloc( sizeof(iterative_solver_c) );
    get_iterative_solver(solver,0);
    solver_init(solver,n,1.0e-10);

    preconditioner_c *pc;
    pc = (preconditioner_c *)malloc( sizeof(preconditioner_c) );
    get_preconditioner(pc,0);
    preconditioner_init(pc,A,0);


    // Fill in some vectors and solve a system
    double dx = 0.01;
    double u[n], f[n];
    for (i=0; i<n; i++) {
        u[i] = 1.0;
        f[i] = 2.0*dx*dx;
    }

    solve(solver,A,u,f,pc,n);


    double ma = 0.0;
    for (i=0; i<n; i++) {
        if (fabs(u[i])>ma) ma = fabs(u[i]);
    }

    if (fabs(ma-0.25)>1.0e-10) {
        printf("-d^2u = 1 should have a maximum of 0.25\n");
        printf("Value found:  %lf\n",ma);
    }*/

    return 0;
}
