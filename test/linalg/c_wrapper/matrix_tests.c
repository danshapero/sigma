#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    void *p;
    int nrow, ncol, mat_type;
} sparse_matrix_c;


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

    // Check to make sure that entries in the matrix are zero
    double z;
    get_value_c(A,0,0,&z);
    printf("%lf \n",z);

    double w=1;
    set_value_c(A,0,0,w);

    get_value_c(A,0,0,&z);
    printf("%lf \n",z);

    return 0;
}
