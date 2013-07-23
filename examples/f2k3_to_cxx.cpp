#include <stdlib.h>
#include <iostream>
#include <string>
#include <time.h>
#include <tetgen.h>
#include <solver.hpp>


int main (int argc, char *argv[]) {

    int i, j, k, n;

    ////////////////
    // Load the mesh
    tetgenio msh;

    char filename[] = "wiffle.1";
    msh.load_tetmesh(filename);

    int nn, nl, ne;
    nn = msh.numberofpoints;
    nl = msh.numberofedges;
    ne = msh.numberoftetrahedra;

    printf("Number of nodes : %d \n",nn);
    printf("Number of edges : %d \n",nl);
    printf("Number of tets  : %d \n",ne);


    //////////////////////////////////
    // Initialize some sparse matrices
    sparse_matrix_c *A;
    A = (sparse_matrix_c *)malloc( sizeof(sparse_matrix_c) );
    get_sparse_matrix_c(A,0);
    init_c(A,nn,nn,nn+2*nl);

    int rows[nn+2*nl], cols[nn+2*nl];
    for (i=0; i<nn; i++) {
        rows[i] = i;
        cols[i] = i;
    }
    for (i=0; i<nl; i++) {
        rows[nn+2*i] = msh.edgelist[2*i]-1;
        cols[nn+2*i] = msh.edgelist[2*i+1]-1;
        rows[nn+2*i+1] = msh.edgelist[2*i+1]-1;
        cols[nn+2*i+1] = msh.edgelist[2*i]-1;
    }
    build_c(A,rows,cols,nn+2*nl);

    int degree = A->max_degree;
    printf("Degree of graph : %d \n",degree);
    int nbrs[degree];
    get_neighbors_c(A,0,nbrs);


    ///////////////////////
    // Fill in the matrices
    int ele[4], d;

    // Iterate through the nodes
    for (i=0; i<nn; i++) {
        d = degree;
        // Find all neighbors of the current node
        get_neighbors_c(A,i,nbrs);

        // Compute the degree of the current node
        for (j=degree-1; j>=0; j--) {
            if (nbrs[j]==-1) d--;
        }

        // Set A[i,i] = degree(i)+1, A[i,j] = -1 for each neighbor j of i
        for (k=0; k<d; k++) {
            j = nbrs[k];
            if (j==i) {
                set_value_c(A,i,i,double(d));
            } else {
                set_value_c(A,i,j,-1.0);
            }
        }
    }

    /*get_neighbors_c(A,0,nbrs,degree);
    double a0[degree];
    d = degree;
    for (j=degree-1; j>=0; j--) {
        if (nbrs[j]==-1) {
            d--;
        } else {
            get_value_c(A,0,nbrs[j],&(a0[j]));
        }
    }
    for (j=0; j<d; j++) printf("%lf ",a0[j]);
    printf("\n");*/


    ////////////////////////////////////
    // Make a solver and preconditioner
    iterative_solver_c *solver;
    solver = (iterative_solver_c *)malloc( sizeof(iterative_solver_c) );
    get_iterative_solver_c(solver,0);
    solver_init_c(solver,nn,1.0e-10);

    preconditioner_c *pc;
    pc = (preconditioner_c *)malloc( sizeof(preconditioner_c) );
    get_preconditioner_c(pc,0);
    preconditioner_init_c(pc,A,0);


    /////////////////////
    // Make some vectors
    double z[2], u[nn], f[nn];
    const double PI = 4.0*atan(1.0);
    srand( time(0) );
    for (i=0; i<nn; i++) {
        z[0] = double(rand())/RAND_MAX;
        z[1] = double(rand())/RAND_MAX;
        f[i] = log(z[0]+1.0e-14)*cos(2*PI*z[1]);
        u[i] = 0.0;
    }


    //////////////////
    // Solve a system
    solve_c(solver,A,u,f,pc,nn);

    return 0;
}


    // tetgenio inp, msh;
    // tetgenbehavior b;
    // char switches[] = "pq1.414a0.1";
    // b.parse_commandline(switches);
    // tetrahedralize(&b, &inp, &msh);
