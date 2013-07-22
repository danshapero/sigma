#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <matrix.h>

int main(int argc, char *argv) {

    // Initialize the matrix
    sparse_matrix_c *A;
    A = (sparse_matrix_c *)malloc( sizeof(sparse_matrix_c) );
    get_sparse_matrix_c(A,0);
    init_c(A,100,100,300);


    // Build the non-zero structure of the matrix
    int rows[300], cols[300];
    int i;
    for (i=0; i<100; i++) {
        rows[i] = i;
        rows[i+100] = i;
        rows[i+200] = i;
        cols[i] = i;
        cols[i+100] = i-1;
        cols[i+200] = i+1;
    }
    cols[100] = 99;
    cols[299] = 0;

    build_c(A,rows,cols,300);


    // Fill in the matrix entries
    for (i=0; i<99; i++) {
        set_value_c(A,i,i,2.0);
        set_value_c(A,i,i+1,-1.0);
        set_value_c(A,i+1,i,-1.0);
    }
    set_value_c(A,99,99,2.0);
    set_value_c(A,0,99,-1.0);
    set_value_c(A,99,0,-1.0);

    double z1, z2, z3;
    get_value_c(A,0,0,&z1);
    get_value_c(A,0,1,&z2);
    get_value_c(A,0,2,&z3);

    if ( !(z1==2.0 && z2==-1.0 && z3==0.0) ) {
        printf("A[0,0:2] should be = [2.0, -1.0, 0.0]\n");
        printf("Values found:  [ %lf, %lf, %lf ]\n",z1,z2,z3);
        return 1;
    }


    // Test matrix multiplication
    double x[100], y[100];
    for (i=0; i<100; i++) {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    matvec_c(A,x,y,100,100);

    double maxval=0.0;
    for (i=0; i<100; i++) {
        if (fabs(y[i])<maxval) {
            maxval = fabs(y[i]);
        }
    }

    if (maxval > 1.0e-14) {
        printf("A*[1.0,...,1.0] should be = 0.0\n");
        printf("Value found:  %lf\n",maxval);
        return 1;
    }

    return 0;
}
