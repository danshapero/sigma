#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char *argv) {

    // Initialize the matrix
    void *A;
    csr_matrix_c(&A);
    csr_init_c(&A,100,100,300);

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

    csr_build_c(&A,rows,cols,300);

    for (i=0; i<99; i++) {
        csr_set_value_c(&A,i,i,2.0);
        csr_set_value_c(&A,i,i+1,-1.0);
        csr_set_value_c(&A,i+1,i,-1.0);
    }
    csr_set_value_c(&A,99,99,2.0);
    csr_set_value_c(&A,0,99,-1.0);
    csr_set_value_c(&A,99,0,-1.0);

    double z1, z2, z3;
    csr_get_value_c(&A,0,0,&z1);
    csr_get_value_c(&A,0,1,&z2);
    csr_get_value_c(&A,0,2,&z3);
    printf("%lf  %lf  %lf\n",z1,z2,z3);

    double x[100], y[100];
    for (i=0; i<100; i++) {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    csr_matvec_c(&A,x,y,100,100);

    double m=0.0;
    for(i=0; i<100; i++) {
        if (fabs(y[i])>m) {
            m = fabs(y[i]);
        }
        //printf("%d   %lf\n",i,y[i]);
    }
    printf("%lf\n",m);

    return 0;
}
