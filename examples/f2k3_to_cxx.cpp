#include <iostream>
#include <string>
#include <tetgen.h>
//#include <matrix.h>
#include <solver.h>

extern"C" {
typedef struct {
    void *p;
    int nrow, ncol, mat_type;
} sparse_matrix_c;

void get_sparse_matrix_c(sparse_matrix_c *cmat, int mat_type);

void init_c(sparse_matrix_c *cmat, int nrow, int ncol, int nnz);

void build_c(sparse_matrix_c *cmat, int *rows, int *cols, int nnz);

void get_value_c(sparse_matrix_c *cmat, int i, int j, double *val);

void set_value_c(sparse_matrix_c *cmat, int i, int j, double val);

void add_value_c(sparse_matrix_c *cmat, int i, int j, double val);

void matvec_c(sparse_matrix_c *cmat, double *x, double *y, int n, int m);
}


int main (int argc, char *argv[]) {

    tetgenio msh;
//    tetgenbehavior b;
//    tetgenio::facet *f;
//    tetgenio::polygon *p;
    int i;

    char filename[] = "wiffle.1";
    msh.load_tetmesh(filename);

//    char switches[] = "pq1.414a0.1";
//    b.parse_commandline(switches);

//    tetrahedralize(&b, &inp, &msh);

    int nn, nl, ne;
    nn = msh.numberofpoints;
    nl = msh.numberofedges;
    ne = msh.numberoftetrahedra;

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

    

    return 0;
}
